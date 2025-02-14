"""
Abstract type representing the evaluation of mass density of each phase (i.e.
units of mass per units of volume, for each cell in the model domain.)
"""
abstract type PhaseMassDensities <: PhaseVariables end


"""
    ConstantCompressibilityDensities(
        sys_or_nph::Union{MultiPhaseSystem, Integer},
        reference_pressure = 1.0,
        reference_density = 0.0,
        compressibility = 1.0
    )

Secondary variable that implements a constant compressibility relationship for
density. Given the reference pressure, compressibility and density at the
reference pressure, each phase density can be computed as:

``ρ(S) = ρ_{ref} e^{(p - p_{ref})c}``

The constructor can take in either one value per phase or a single value for all
phases for the reference pressure, compressibility and density at reference
conditions.

## Fields
$FIELDS
"""
struct ConstantCompressibilityDensities{T} <: PhaseMassDensities
    "Reference pressure for each phase (where the reference densities are given)"
    reference_pressure::T
    "Densities at the reference point"
    reference_densities::T
    "Compressibility factor used when expanding around reference pressure, typically between 1e-3 and 1e-10"
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

function ConstantCompressibilityDensities(; p_ref = DEFAULT_MINIMUM_PRESSURE, density_ref = 1000.0, compressibility = 1e-10)
    n = max(length(p_ref), length(density_ref), length(compressibility))
    return ConstantCompressibilityDensities(n, p_ref, density_ref, compressibility)
end

@jutul_secondary function update_density!(rho, density::ConstantCompressibilityDensities, model, Pressure, ix)
    p_ref, c, rho_ref = density.reference_pressure, density.compressibility, density.reference_densities
    for i in ix
        for ph in axes(rho, 1)
            @inbounds rho[ph, i] = constant_expansion(Pressure[i], p_ref[ph], c[ph], rho_ref[ph])
        end
    end
end

@inline function constant_expansion(p::Real, p_ref::Real, c::Real, f_ref::Real)
    Δ = p - p_ref
    return f_ref * exp(Δ * c)
end

# Total masses
@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::SimulationModel{G, S}, PhaseMassDensities, FluidVolume, ix) where {G, S<:SinglePhaseSystem}
    @inbounds for i in ix
        V = FluidVolume[i]
        @inbounds for ph in axes(totmass, 1)
            totmass[ph, i] = PhaseMassDensities[ph, i]*V
        end
    end
end

@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::SimulationModel{G, S}, PhaseMassDensities, Saturations, FluidVolume, ix) where {G, S<:ImmiscibleSystem}
    rho = PhaseMassDensities
    s = Saturations
    @inbounds for i in ix
        V = FluidVolume[i]
        @inbounds for ph in axes(totmass, 1)
            totmass[ph, i] = rho[ph, i]*V*s[ph, i]
        end
    end
end

# Total mass
@jutul_secondary function update_total_mass!(totmass, tv::TotalMass, model::SimulationModel{G, S}, TotalMasses, ix) where {G, S<:MultiPhaseSystem}
    @inbounds for c in ix
        tmp = zero(eltype(totmass))
        @inbounds for ph in axes(TotalMasses, 1)
            tmp += TotalMasses[ph, i]
        end
        totmass[c] = tmp
    end
end

@jutul_secondary function update_phase_mass_mob!(ρλ, var::PhaseMassMobilities, model, PhaseMassDensities, PhaseMobilities, ix)
    for i in ix
        @inbounds for ph in axes(ρλ, 1)
            ρλ[ph, i] = PhaseMassDensities[ph, i]*PhaseMobilities[ph, i]
        end
    end
end


@jutul_secondary function update_phase_mass_mob!(λ, var::PhaseMobilities, model, RelativePermeabilities, PhaseViscosities, ix)
    for i in ix
        @inbounds for ph in axes(λ, 1)
            λ[ph, i] = RelativePermeabilities[ph, i]/PhaseViscosities[ph, i]
        end
    end
end

"""
    FluidVolume()

Variable typically taken to be a parameter. Represents the per-cell volume that
where multiphase flow can occur. For a well, this is the volume inside the
well-bore free flow can occur. For a porous medium, this is the void space
inside the pores that is, to some extent, connected and open to flow (effective
pore-volume).
"""
struct FluidVolume <: ScalarVariable end
Jutul.minimum_value(::FluidVolume) = eps()

function Jutul.default_parameter_values(data_domain, model, param::FluidVolume, symb)
    vol = missing
    if haskey(data_domain, :fluid_volume, Cells())
        vol = data_domain[:fluid_volume]
    else
        vol = pore_volume(data_domain, throw = false)
    end
    if ismissing(vol)
        g = physical_representation(data_domain)
        vol = domain_fluid_volume(g)
        if ismissing(vol)
            error(":volumes or :pore_volume symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
        end
    end
    return copy(vol)
end

Base.@kwdef struct Temperature{T} <: ScalarVariable
    min::T = 273.15
    max::T = 1e6
    max_rel::Union{T, Nothing} = nothing
    max_abs::Union{T, Nothing} = nothing
end

Jutul.default_value(model, T::Temperature) = 303.15 # 30.15 C°
function Jutul.default_parameter_values(data_domain, model, param::Temperature, symb)
    T_default = Jutul.default_value(model, param)
    nc = number_of_cells(data_domain)
    if haskey(data_domain, :temperature)
        T_domain = data_domain[:temperature]
        @assert length(T_domain) == nc
        T = copy(T_domain)
    else
        T = fill(T_default, nc)
    end
    return T
end

Jutul.minimum_value(T::Temperature) = T.min
Jutul.maximum_value(T::Temperature) = T.max
Jutul.absolute_increment_limit(T::Temperature) = T.max_abs
Jutul.relative_increment_limit(T::Temperature) = T.max_rel
