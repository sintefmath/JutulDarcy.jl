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

struct BrineCO2MixingDensities{T} <: AbstractCompositionalDensities
    tab::T
    coeffs::NTuple{4, Float64}
end

function BrineCO2MixingDensities(tab::T; coefficients = (37.51, −9.585*10e−2, 8.74*10e−4, −5.044*10e−7)) where T
    return BrineCO2MixingDensities{T}(tab, coefficients)
end

@jutul_secondary function update_density!(rho, rho_def::BrineCO2MixingDensities, model::SimulationModel{D, S}, Pressure, Temperature, LiquidMassFractions, ix) where {D, S<:CompositionalSystem}
    # V = 10e−6*(37.51 − 9.585*10e−2*T + 8.74*10e−4*T^2 − 5.044*10e−7*T^3)
    # M_co2 = 44.01*10^-3
    # ρ_co2 = M_co2/V

    c1, c2, c3, c4 = rho_def.coeffs
    sys = model.system
    eos = sys.equation_of_state
    cnames = eos.mixture.component_names
    # TODO: This is hard coded.
    @assert cnames[1] == raw"CarbonDioxide"
    @assert length(cnames) == 2
    l, v = phase_indices(sys)
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        X_co2 = LiquidMassFractions[1, i]
        X_h2o = LiquidMassFractions[2, i]

        rho_co2, rho_h2o_pure = rho_def.tab(p, T)

        vol_co2 = 10e−6*(c1 + c2*T + c3*T^2 + c4*T^3)
        rho_liquid_co2_pure = 44.01*10^-3/vol_co2
        rho_brine_mix = rho_liquid_co2_pure*X_co2 + rho_h2o_pure*X_h2o

        rho[v, i] = rho_co2
        rho[l, i] = rho_brine_mix
    end
end
