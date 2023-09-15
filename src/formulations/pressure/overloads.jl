function Jutul.select_secondary_variables!(pvar, ::PressureFormulation, model::PressureModel)
    pvar[:PressureReductionFactors] = PressureReductionFactors()
end

function Jutul.select_equations!(eqs, ::PressureFormulation, model::PressureModel)
    @assert haskey(eqs, :mass_conservation)
    eqs[:pressure] = PressureEquation(eqs[:mass_conservation])
    delete!(eqs, :mass_conservation)
end

function Jutul.setup_equation_storage(model, eq::PressureEquation, storage; extra_sparsity = nothing, kwarg...)
    s = ConservationLawTPFAStorage(model, eq.conservation; kwarg...)
    w = storage[:state][:PressureReductionFactors]
    return (conservation = s, weights = w)
end

function Jutul.update_equation!(eq_s, p::PressureEquation, storage, model, dt)
    Jutul.update_equation!(eq_s.conservation, p.conservation, storage, model, dt)
end

function Jutul.declare_pattern(model, p::PressureEquation, e_s, entity)
    Jutul.declare_pattern(model, p.conservation, e_s.conservation, entity)
end

function Jutul.align_to_jacobian!(p_eq_s, p_eq::PressureEquation, jac, model, u::Cells; equation_offset = 0, variable_offset = 0)
    eq_s = p_eq_s.conservation
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
function Jutul.get_diagonal_entries(eq_p::PressureEquation, eq_s)
    return Jutul.get_diagonal_entries(eq_p.conservation, eq_s.conservation)
end

function Jutul.update_linearized_system_equation!(nz, r, model, p::PressureEquation, peq_s)
    eq_s = peq_s.conservation
    w = peq_s.weights
    law = p.conservation
    acc = eq_s.accumulation
    cell_flux = eq_s.half_face_flux_cells
    face_flux = eq_s.half_face_flux_faces
    conn_pos = law.flow_discretization.conn_pos
    nc, ncomp, np = Jutul.ad_dims(acc)
    @assert np == 1
    @assert !Jutul.use_sparse_sources(law)
    @assert isnothing(face_flux)

    Ncomp = Val(ncomp)
    for cell = 1:nc
        fill_pressure_eq_cell!(nz, r, cell, w, acc, cell_flux, conn_pos, Ncomp)
    end
end

function fill_pressure_eq_cell!(nz, r, cell, weights::Matrix{T}, acc, cell_flux, conn_pos, ::Val{N_comp}) where {N_comp, T}
    eq = zero(T)
    # Sum up accumulation w_i*∂M_i
    for e in 1:N_comp
        a_e = Jutul.get_entry(acc, cell, e)
        eq += a_e*weights[e, cell]
    end
    # Sum up the fluxes
    start = conn_pos[cell]
    stop = conn_pos[cell + 1] - 1
    for i = start:stop
        v_p = zero(T)
        for e in 1:N_comp
            v_p += weights[e, cell]*Jutul.get_entry(cell_flux, i, e)
        end
        eq -= v_p

        fpos = Jutul.get_jacobian_pos(cell_flux, i, 1, 1)
        if fpos > 0
            Jutul.update_jacobian_inner!(nz, fpos, only(v_p.partials))
        end
    end
    apos = Jutul.get_jacobian_pos(acc, cell, 1, 1)
    Jutul.update_jacobian_inner!(nz, apos, only(eq.partials))
    r[cell] = eq.value
end

function Jutul.apply_forces_to_equation!(acc, storage, model::PressureModel, eq::PressureEquation, eq_s, force, time)
    Jutul.apply_forces_to_equation!(acc, storage, model, eq.conservation, eq_s.conservation, force, time)
end

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
