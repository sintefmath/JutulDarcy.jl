function Jutul.convergence_criterion(model::SimpleWellFlowModel, storage, eq::ConservationLaw{:TotalMasses}, eq_s, r; dt = 1)
    vol = value(only(storage.state.FluidVolume))
    scale = 0.1
    e = map(x -> scale*abs(x)*dt/vol, vec(r))
    R = (CNV = (errors = e, names = map(x -> "M$x", eachindex(e))), )
    return R
end

@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::SimpleWellFlowModel,
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
    @. eq_buf = (M - M₀)/Δt
end
