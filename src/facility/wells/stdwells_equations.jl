function Jutul.setup_storage_system!(storage, model::StandardWellFlowModel, system::SimpleWellSystem)
    number_of_residuals = length(storage.state.TotalMasses)
    storage[:well_convergence_names] = ["M$i" for i in 1:number_of_residuals]
    return storage
end

function Jutul.convergence_criterion(model::StandardWellFlowModel, storage, eq::ConservationLaw{:TotalMasses}, eq_s, r; dt = 1, update_report = missing)
    vol = storage.state.FluidVolume
    scale = 0.1
    # Reports retain the errors from each iteration, so keep one copy of the
    # residual and perform all scaling on that copy.
    e = copy(vec(r))
    if ndims(r) == 1
        @. e = scale*abs(e)*dt/value(vol)
    else
        errors_by_cell = reshape(e, size(r))
        volumes_by_cell = reshape(vol, 1, :)
        @. errors_by_cell = scale*abs(errors_by_cell)*dt/value(volumes_by_cell)
    end
    R = (CNV = (errors = e, names = storage.well_convergence_names), )
    return R
end

@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::StandardWellFlowModel,
    Pressure,
    MassFractions,
    FluidVolume,
    ix)
    density(P) = 1.0 + (P - DEFAULT_MINIMUM_PRESSURE)*1e-8

    for cell in ix
        rho = density(Pressure[cell])
        vol = FluidVolume[cell]
        V = vol*rho
        for i in axes(MassFractions, 1)
            totmass[i, cell] = V*MassFractions[i, cell]
        end
    end
end

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::ConservationLaw, model::SimpleWellFlowModel, Δt, ldisc = local_discretization(eq, self_cell)) where T_e
    conserved = conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]
    if ndims(M) == 1
        eq_buf[1] = (M[self_cell] - M₀[self_cell])/Δt
    else
        for i in eachindex(eq_buf)
            eq_buf[i] = (M[i, self_cell] - M₀[i, self_cell])/Δt
        end
    end
    return eq_buf
end

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::ConservationLaw, model::SimpleWellModel, Δt, ldisc = local_discretization(eq, self_cell)) where T_e
    conserved = conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]
    if ndims(M) == 1
        eq_buf[1] = (M[self_cell] - M₀[self_cell])/Δt
    else
        for i in eachindex(eq_buf)
            eq_buf[i] = (M[i, self_cell] - M₀[i, self_cell])/Δt
        end
    end
    return eq_buf
end
