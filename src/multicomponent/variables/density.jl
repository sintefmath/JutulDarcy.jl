abstract type AbstractCompositionalDensities <: PhaseMassDensities end

struct TwoPhaseCompositionalDensities <: AbstractCompositionalDensities end

struct ThreePhaseCompositionalDensitiesLV{T} <: AbstractCompositionalDensities
    immiscible_pvt::T
end


@jutul_secondary function update_density!(rho, m::TwoPhaseCompositionalDensities, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ix) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    l, v = phase_indices(sys)
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        rho[l, i], rho[v, i] = mass_densities(eos, p, T, FlashResults[i])
    end
end

@jutul_secondary function update_density!(rho, m::ThreePhaseCompositionalDensitiesLV, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ix) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    rhos = reference_densities(sys)

    pvt = m.immiscible_pvt
    a, l, v = phase_indices(sys)
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        rho[l, i], rho[v, i] = mass_densities(eos, p, T, FlashResults[i])
        rho[a, i] = rhos[a]*shrinkage(pvt, p)
    end
end
