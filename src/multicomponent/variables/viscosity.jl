abstract type AbstractCompositionalViscosities <: PhaseVariables end

struct LBCViscosities <: AbstractCompositionalViscosities end

struct ThreePhaseLBCViscositiesLV{T} <: AbstractCompositionalViscosities
    immiscible_pvt::T
end


@jutul_secondary function update_as_secondary!(mu, m::LBCViscosities, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    n = size(mu, 2)
    tb = minbatch(model.context)
    l, v = phase_indices(sys)
    @inbounds @batch minbatch = tb for i in 1:n
        p = Pressure[i]
        T = Temperature[i]
        mu[l, i], mu[v, i] = lbc_viscosities(eos, p, T, FlashResults[i])
    end
end

@jutul_secondary function update_as_secondary!(mu, m::ThreePhaseLBCViscositiesLV, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    pvt = m.immiscible_pvt
    n = size(mu, 2)
    tb = minbatch(model.context)
    a, l, v = phase_indices(sys)
    @inbounds @batch minbatch = tb for i in 1:n
        p = Pressure[i]
        T = Temperature[i]
        mu[l, i], mu[v, i] = lbc_viscosities(eos, p, T, FlashResults[i])
        mu[a, i] = viscosity(pvt, p)
    end
end
