abstract type AbstractCompositionalViscosities <: PhaseVariables end

struct LBCViscosities <: AbstractCompositionalViscosities end

struct ThreePhaseLBCViscositiesLV{T} <: AbstractCompositionalViscosities
    immiscible_pvt::T
end


@jutul_secondary function update_viscosity!(mu, m::LBCViscosities, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ix) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    l, v = phase_indices(sys)
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        mu[l, i], mu[v, i] = lbc_viscosities(eos, p, T, FlashResults[i])
    end
end

@jutul_secondary function update_viscosity!(mu, m::ThreePhaseLBCViscositiesLV, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ix) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    pvt = m.immiscible_pvt
    a, l, v = phase_indices(sys)
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        mu[l, i], mu[v, i] = lbc_viscosities(eos, p, T, FlashResults[i])
        mu[a, i] = viscosity(pvt, p)
    end
end

struct PTViscosities{T} <: AbstractCompositionalViscosities
    tab::T
end

@jutul_secondary function update_viscosity!(mu, mu_def::PTViscosities, model, Pressure, Temperature, ix)
    tab = mu_def.tab
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        μ = tab(p, T)
        @. mu[:, i] = μ
    end
end
