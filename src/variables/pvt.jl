"""
Mass density of each phase
"""
abstract type PhaseMassDensities <: PhaseVariables end

struct ConstantCompressibilityDensities{T} <: PhaseMassDensities
    reference_pressure::T
    reference_densities::T
    compressibility::T
    function ConstantCompressibilityDensities(sys_or_nph::Union{MultiPhaseSystem, Integer}, reference_pressure = 101325.0, reference_density = 1000.0, compressibility = 1e-10)
        if isa(sys_or_nph, Integer)
            nph = sys_or_nph
        else
            nph = number_of_phases(sys_or_nph)
        end

        pref = expand_to_phases(reference_pressure, nph)
        rhoref = expand_to_phases(reference_density, nph)
        c = expand_to_phases(compressibility, nph)
        T = typeof(c)
        new{T}(pref, rhoref, c)
    end
end

function Base.show(io::IO, t::MIME"text/plain", d::ConstantCompressibilityDensities)
    p_r = d.reference_pressure./1e5
    ρ_r = d.reference_densities
    print(io, "ConstantCompressibilityDensities (ref_dens=$ρ_r kg/m^3, ref_p=$p_r bar)")
end

function ConstantCompressibilityDensities(; p_ref = 101325.0, density_ref = 1000.0, compressibility = 1e-10)
    n = max(length(p_ref), length(density_ref), length(compressibility))
    return ConstantCompressibilityDensities(n, p_ref, density_ref, compressibility)
end

@jutul_secondary function update_as_secondary!(rho, density::ConstantCompressibilityDensities, model, Pressure)
    p_ref, c, rho_ref = density.reference_pressure, density.compressibility, density.reference_densities
    @tullio rho[ph, i] = constant_expansion(Pressure[i], p_ref[ph], c[ph], rho_ref[ph])
end

@inline function constant_expansion(p::Real, p_ref::Real, c::Real, f_ref::Real)
    Δ = p - p_ref
    return f_ref * exp(Δ * c)
end

# Total masses
@jutul_secondary function update_as_secondary!(totmass, tv::TotalMasses, model::SimulationModel{G, S}, PhaseMassDensities, FluidVolume) where {G, S<:SinglePhaseSystem}
    @tullio totmass[ph, i] = PhaseMassDensities[ph, i]*FluidVolume[i]
end

@jutul_secondary function update_as_secondary!(totmass, tv::TotalMasses, model::SimulationModel{G, S}, PhaseMassDensities, Saturations, FluidVolume) where {G, S<:ImmiscibleSystem}
    rho = PhaseMassDensities
    s = Saturations
    @tullio totmass[ph, i] = rho[ph, i]*FluidVolume[i]*s[ph, i]
end

# Total mass
@jutul_secondary function update_as_secondary!(totmass, tv::TotalMass, model::SimulationModel{G, S}, TotalMasses) where {G, S<:MultiPhaseSystem}
    @tullio totmass[i] = TotalMasses[ph, i]
end

@jutul_secondary function update_as_secondary!(ρλ, var::PhaseMassMobilities, model, RelativePermeabilities, PhaseMassDensities, PhaseViscosities)
    mb = minbatch(model.context)
    @batch minbatch = mb for i in axes(ρλ, 2)
        @inbounds for ph in axes(ρλ, 1)
            ρλ[ph, i] = PhaseMassDensities[ph, i]*RelativePermeabilities[ph, i]/PhaseViscosities[ph, i]
        end
    end
end

struct FluidVolume <: ScalarVariable end
Jutul.default_values(model, ::FluidVolume) = fluid_volume(model.domain)
Jutul.minimum_value(::FluidVolume) = eps()

struct Temperature <: ScalarVariable end

Jutul.default_value(model, ::Temperature) = 303.15 # 30.15 C°
