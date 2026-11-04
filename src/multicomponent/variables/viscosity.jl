abstract type AbstractCompositionalViscosities <: PhaseVariables end

struct LBCViscosities{C} <: AbstractCompositionalViscosities
    coefficients::C
end

LBCViscosities(; coefficients = nothing) = LBCViscosities(coefficients)

struct ThreePhaseLBCViscositiesLV{T, C} <: AbstractCompositionalViscosities
    immiscible_pvt::T
    coefficients::C
end

ThreePhaseLBCViscositiesLV(pvt; coefficients = nothing) = ThreePhaseLBCViscositiesLV(pvt, coefficients)

compositional_lbc_viscosities(eos, p, T, flash, ::Nothing) = lbc_viscosities(eos, p, T, flash)

function compositional_lbc_viscosities(eos, p, T, flash, coefficients)
    has_liquid = MultiComponentFlash.liquid_phase_present(flash.state)
    has_vapor = MultiComponentFlash.vapor_phase_present(flash.state)
    liquid = has_liquid ? flash.liquid : flash.vapor
    vapor = has_vapor ? flash.vapor : flash.liquid
    mu_l = lbc_viscosity(eos, p, T, liquid; coeff = coefficients)
    mu_v = lbc_viscosity(eos, p, T, vapor; coeff = coefficients)
    return (μ_l = mu_l, μ_v = mu_v)
end

@jutul_secondary function update_viscosity!(mu, m::LBCViscosities, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ix) where {D, S<:CompositionalSystem}
    sys = model.system
    eos = sys.equation_of_state
    l, v = phase_indices(sys)
    @inbounds for i in ix
        p = Pressure[i]
        T = Temperature[i]
        mu[l, i], mu[v, i] = compositional_lbc_viscosities(eos, p, T, FlashResults[i], m.coefficients)
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
        mu[l, i], mu[v, i] = compositional_lbc_viscosities(eos, p, T, FlashResults[i], m.coefficients)
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
        for ph in eachindex(μ)
            mu[ph, i] = μ[ph]
        end
    end
end
