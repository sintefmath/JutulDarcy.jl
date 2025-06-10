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

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::PressureEquation, model, Δt, ldisc = local_discretization(eq, self_cell)) where T_e
    # Compute accumulation term
    ceq = eq.conservation
    conserved = conserved_symbol(ceq)
    M₀ = state0[conserved]
    M = state[conserved]
    w = state.PressureReductionFactors
    # Compute ∇⋅V
    if false
        disc = ceq.flow_discretization
        flux(face) = Jutul.face_flux(face, ceq, state, model, Δt, disc, ldisc, Val(T_e))
        div_v = ldisc.div(flux)
        val = zero(eltype(eq_buf))
        for i in eachindex(div_v)
            ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, i, self_cell)
            w_i = w[i, self_cell]
            @inbounds val += w_i*(∂M∂t + div_v[i])
        end
    else
        tmp = zeros(eltype(eq_buf), size(M, 1))
        Jutul.update_equation_in_entity!(tmp, self_cell, state, state0, ceq, model, Δt, ldisc)
        val = sum(w[:, self_cell].*tmp)
    end
    @assert length(eq_buf) == 1
    # error("Generic AD pressure equation is not yet functional")
    eq_buf[1] = val
end

function Jutul.update_equation!(eq_s::PressureEquationTPFAStorage, eq_p::PressureEquation, storage, model, dt)
    pressure_update_accumulation!(eq_s, eq_p, storage, model, dt)
    pressure_update_half_face_flux!(eq_s, eq_p, storage.state, model, dt, eq_p.conservation.flow_discretization)
end

function pressure_update_accumulation!(eq_s, eq_p, storage, model, dt)
    conserved = eq_s.accumulation_symbol
    acc = Jutul.get_entries(eq_s.accumulation)
    m0, m = Jutul.state_pair(storage, conserved, model)
    w = storage.state.PressureReductionFactors
    pressure_update_accumulation_inner!(acc, m, m0, dt, w)
    return acc
end

function pressure_update_accumulation_inner!(acc, m, m0, dt, w)
    for cell in axes(m, 2)
        val = zero(eltype(acc))
        for i in axes(m, 1)
            val += w[i, cell]*Jutul.accumulation_term(m, m0, dt, i, cell)
        end
        acc[1, cell] = val
    end
    return acc
end

function pressure_update_half_face_flux!(eq_s::PressureEquationTPFAStorage, eq_p::PressureEquation, state, model, dt, flow_disc)
    flux_c = get_entries(eq_s.half_face_flux_cells)

    T = eltype(flux_c)
    N = number_of_components(model.system)
    zero_flux = zero(SVector{N, T})
    state_c = Jutul.local_ad(state, 1, T)
    w = state.PressureReductionFactors
    pressure_update_half_face_flux_tpfa!(flux_c, zero_flux, eq_p, state_c, w, model, dt, flow_disc, Cells())
end

function pressure_update_half_face_flux_tpfa!(hf_cells, zero_flux::SVector{N, T}, eq, state::S, w, model, dt, flow_disc, ::Cells) where {T, N, S<:LocalStateAD}
    conn_data = flow_disc.conn_data
    conn_pos = flow_disc.conn_pos
    M = global_map(model.domain)
    nc = length(conn_pos)-1
    eq_cons = eq.conservation
    @tic "flux (cells)" for c in 1:nc
        self = Jutul.full_cell(c, M)
        state_c = Jutul.new_entity_index(state, self)
        pressure_update_half_face_flux_tpfa_internal!(hf_cells, zero_flux, eq_cons, state_c, w, model, dt, flow_disc, conn_pos, conn_data, c)
    end
end

function pressure_update_half_face_flux_tpfa_internal!(hf_cells, zero_flux, eq, state, w, model, dt, flow_disc, conn_pos, conn_data, c)
    start = @inbounds conn_pos[c]
    stop = @inbounds conn_pos[c+1]-1
    for i in start:stop
        (; self, other, face, face_sign) = @inbounds conn_data[i]
        @assert self == c
        q = Jutul.face_flux!(zero_flux, self, other, face, face_sign, eq, state, model, dt, flow_disc)
        val = zero(eltype(hf_cells))
        for j in 1:length(q)
            val += q[j]*w[j, c]
        end
        @inbounds hf_cells[i] = val
    end
end

function Jutul.update_linearized_system_equation!(nz, r, model, peq::PressureEquation, eq_s::PressureEquationTPFAStorage)
    acc = eq_s.accumulation
    cell_flux = eq_s.half_face_flux_cells
    cpos = peq.conservation.flow_discretization.conn_pos
    ctx = model.context
    Jutul.update_linearized_system_subset_conservation_accumulation!(nz, r, model, acc, cell_flux, cpos, ctx)
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

function Jutul.get_diagonal_entries(eq::PressureEquation, eq_s::PressureEquationTPFAStorage)
    return eq_s.accumulation.entries
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

function Jutul.convergence_criterion(model, storage, eq::PressureEquation, eq_s, r; dt = 1.0, update_report = missing)
    # M = global_map(model.domain)
    # v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    # Φ = v(storage.state.FluidVolume)
    # p_res = dt*sum(abs, r)/sum(Φ)
    # R = (Residual = (errors = (p_res,), names = (:L1, )),)

    dp_abs, dp_rel = JutulDarcy.pressure_increments(model, storage.state, update_report)
    R = (
        increment_dp_abs = (errors = (dp_abs/1e6, ), names = (raw"Δp (abs, MPa)", ), ),
        increment_dp_rel = (errors = (dp_rel, ), names = (raw"Δp (rel)", ), ),
    )
    return R
end

@jutul_secondary function update_rs!(rs, ph::Rs, model::SimulationModel{D, S, F}, Pressure, BlackOilUnknown, ix)  where {D, S<:BlackOilVariableSwitchingSystem, F<:PressureFormulation}
    T = eltype(rs)
    @inbounds for i in ix
        X = BlackOilUnknown[i]
        phases = X.phases_present
        p = @inbounds Pressure[i]
        reg_i = JutulDarcy.region(ph.regions, i)
        rs_max = JutulDarcy.table_by_region(model.system.rs_max, reg_i)
        rs_max_val = rs_max(p)

        if phases == JutulDarcy.OilOnly
            r = X.val
            if r >= rs_max_val
                r = rs_max_val
            end
        else
            r = rs_max_val
        end
        rs[i] = r
    end
end
