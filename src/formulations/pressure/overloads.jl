function Jutul.select_secondary_variables!(pvar, ::PressureFormulation, model::PressureModel)
    pvar[:PressureReductionFactors] = PressureReductionFactors()
end

function Jutul.select_equations!(eqs, ::PressureFormulation, model::PressureModel)
    @assert haskey(eqs, :mass_conservation)
    eqs[:pressure] = PressureEquation(eqs[:mass_conservation])
    delete!(eqs, :mass_conservation)
end

function Jutul.local_discretization(eq::PressureEquation, self_cell)
    return Jutul.local_discretization(eq.conservation, self_cell)
end

function update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::PressureEquation, model, Δt, ldisc = local_discretization(eq, self_cell)) where T_e
    # Compute accumulation term
    ceq = eq.conservation
    conserved = conserved_symbol(ceq)
    M₀ = state0[conserved]
    M = state[conserved]
    w = state.PressureReductionFactors
    # Compute ∇⋅V
    disc = ceq.flow_discretization
    flux(face) = Jutul.face_flux(face, ceq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)
    eq_buf[1] = zero(eltype(eq_buf))
    for i in eachindex(div_v)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, i, self_cell)
        w_i = w[i, self_cell]
        @inbounds eq_buf[1] += w_i*(∂M∂t + div_v[i])
    end
end

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::PressureEquation, model::SimpleWellModel, Δt, ldisc = local_discretization(eq, self_cell)) where T_e
    w = state.PressureReductionFactors
    @assert size(w, 2) == 1
    conserved = conserved_symbol(eq.conservation)
    M₀ = state0[conserved]
    M = state[conserved]
    val = zero(eltype(eq_buf))
    for i in axes(M, 1)
        val += w[i, 1]*(M[i, 1] - M₀[i, 1])
    end
    eq_buf[1] = val/Δt
end

struct PressureEquationTPFAStorage{A, HC}
    accumulation::A
    accumulation_symbol::Symbol
    half_face_flux_cells::HC
end

function PressureEquationTPFAStorage(model, eq::PressureEquation; kwarg...)
    ceq = eq.conservation
    ceq.flow_discretization::TwoPointPotentialFlowHardCoded
    number_of_equations = 1
    D, ctx = model.domain, model.context
    cell_entity = Cells()
    face_entity = Faces()
    nc = count_active_entities(D, cell_entity, for_variables = false)
    nf = count_active_entities(D, face_entity, for_variables = false)
    nhf = number_of_half_faces(ceq.flow_discretization)
    face_partials = degrees_of_freedom_per_entity(model, face_entity)
    @assert face_partials == 0 "Only supported for cell-centered discretization"
    alloc = (n, entity, n_entities_pos) -> CompactAutoDiffCache(number_of_equations, n, model,
                                                                                entity = entity, n_entities_pos = n_entities_pos, 
                                                                                context = ctx; kwarg...)
    # Accumulation terms
    acc = alloc(nc, cell_entity, nc)
    # Source terms - as sparse matrix
    t_acc = eltype(acc.entries)
    # src = sparse(zeros(0), zeros(0), zeros(t_acc, 0), size(acc.entries)...)
    # Half face fluxes - differentiated with respect to pairs of cells
    hf_cells = alloc(nhf, cell_entity, nhf)
    # # Half face fluxes - differentiated with respect to the faces
    # if face_partials > 0
    #     hf_faces = alloc(nf, face_entity, nhf)
    # else
    #     hf_faces = nothing
    # end
    return PressureEquationTPFAStorage(acc, conserved_symbol(ceq), hf_cells)
end

function Jutul.setup_equation_storage(model,
        eq::PressureEquation{ConservationLaw{A, B, C, D}}, storage; extra_sparsity = nothing, kwarg...
        ) where {A, B<:TwoPointPotentialFlowHardCoded, C, D}
    return PressureEquationTPFAStorage(model, eq; kwarg...)
end

# function Jutul.update_equation!(eq_s, p::PressureEquation, storage, model, dt)
#     Jutul.update_equation!(eq_s.conservation, p.conservation, storage, model, dt)
# end

function Jutul.declare_pattern(model, peq::PressureEquation{ConservationLaw{A, B, C, D}}, e_s, entity::Cells) where {A, B<:TwoPointPotentialFlowHardCoded, C, D}
    df = peq.conservation.flow_discretization
    hfd = Array(df.conn_data)
    n = number_of_entities(model, peq)
    # Diagonals
    diagonals = [i for i in 1:n]
    if length(hfd) > 0
        # Fluxes
        I = map(x -> x.self, hfd)
        J = map(x -> x.other, hfd)
        I, J = Jutul.map_ij_to_active(I, J, model.domain, entity)
        I = vcat(I, diagonals)
        J = vcat(J, diagonals)
    else
        I = collect(diagonals)
        J = collect(diagonals)
    end
    return (I, J)
end

function Jutul.align_to_jacobian!(eq_s::PressureEquationTPFAStorage, p_eq::PressureEquation, jac, model, u::Cells; equation_offset = 0, variable_offset = 0)
    eq = p_eq.conservation
    fd = eq.flow_discretization
    M = global_map(model.domain)

    acc = eq_s.accumulation
    hflux_cells = eq_s.half_face_flux_cells
    nu, = Jutul.ad_dims(acc)
    dims = (nu, 1, 1) # Assume 1 eq, 1 pressure atm
    Jutul.diagonal_alignment!(acc, eq, jac, u, model.context,
        target_offset = equation_offset,
        source_offset = variable_offset,
        dims = dims
    )
    Jutul.half_face_flux_cells_alignment!(hflux_cells, acc, jac, model.context, M, fd,
        target_offset = equation_offset,
        source_offset = variable_offset,
        dims = dims
    )
end

# function Jutul.get_diagonal_entries(eq_p::PressureEquation, eq_s)
#     return Jutul.get_diagonal_entries(eq_p.conservation, eq_s.conservation)
# end

# function Jutul.update_linearized_system_equation!(nz, r, model, p::PressureEquation, peq_s)
#     eq_s = peq_s.conservation
#     w = peq_s.weights
#     law = p.conservation
#     acc = eq_s.accumulation
#     cell_flux = eq_s.half_face_flux_cells
#     face_flux = eq_s.half_face_flux_faces
#     conn_pos = law.flow_discretization.conn_pos
#     nc, ncomp, np = Jutul.ad_dims(acc)
#     @assert np == 1
#     @assert !Jutul.use_sparse_sources(law)
#     @assert isnothing(face_flux)

#     Ncomp = Val(ncomp)
#     for cell = 1:nc
#         fill_pressure_eq_cell!(nz, r, cell, w, acc, cell_flux, conn_pos, Ncomp)
#     end
# end

# function fill_pressure_eq_cell!(nz, r, cell, weights::Matrix{T}, acc, cell_flux, conn_pos, ::Val{N_comp}) where {N_comp, T}
#     eq = zero(T)
#     # Sum up accumulation w_i*∂M_i
#     for e in 1:N_comp
#         a_e = Jutul.get_entry(acc, cell, e)
#         eq += a_e*weights[e, cell]
#     end
#     # Sum up the fluxes
#     start = conn_pos[cell]
#     stop = conn_pos[cell + 1] - 1
#     for i = start:stop
#         v_p = zero(T)
#         for e in 1:N_comp
#             v_p -= weights[e, cell]*Jutul.get_entry(cell_flux, i, e)
#         end
#         eq -= v_p

#         fpos = Jutul.get_jacobian_pos(cell_flux, i, 1, 1)
#         if fpos > 0
#             Jutul.update_jacobian_inner!(nz, fpos, only(v_p.partials))
#         end
#     end
#     apos = Jutul.get_jacobian_pos(acc, cell, 1, 1)
#     Jutul.update_jacobian_inner!(nz, apos, only(eq.partials))
#     r[cell] = eq.value
# end

# function Jutul.apply_forces_to_equation!(acc, storage, model::PressureModel, eq::PressureEquation, eq_s, force, time)
#     Jutul.apply_forces_to_equation!(acc, storage, model, eq.conservation, eq_s.conservation, force, time)
# end

function Jutul.select_primary_variables!(pvar, ::PressureFormulation, model::PressureModel)
    for k in keys(pvar)
        if k != :Pressure
            delete!(pvar, k)
        end
    end
end

function convergence_criterion(model, storage, eq::PressureEquation, eq_s, r; dt = 1.0, update_report = missing)
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = v(storage.state.FluidVolume)

    p_res = dt*sum(abs, r)/sum(Φ)
    R = (Residual = (errors = (p_res,), names = (:L1, )),)
    return R
end
