export PhaseMassDensities, ConstantCompressibilityDensities
export BrooksCoreyRelPerm, TabulatedRelPermSimple

abstract type AbstractRelativePermeabilities <: PhaseVariables end
struct RelativePermeabilities <: AbstractRelativePermeabilities end
struct FluidVolume <: ScalarVariable end
Jutul.default_values(model, ::FluidVolume) = fluid_volume(model.domain)

struct PhaseViscosities <: PhaseVariables end
Jutul.default_value(model, v::PhaseViscosities) = 1e-3

degrees_of_freedom_per_entity(model, sf::PhaseVariables) = number_of_phases(model.system)

# Single-phase specialization
degrees_of_freedom_per_entity(model::SimulationModel{D, S}, sf::ComponentVariable) where {D, S<:SinglePhaseSystem} = 1

# Immiscible specialization
degrees_of_freedom_per_entity(model::SimulationModel{D, S}, sf::ComponentVariable) where {D, S<:ImmiscibleSystem} = number_of_phases(model.system)

function select_secondary_variables!(S, system::MultiPhaseSystem, model)
    select_default_darcy_secondary_variables!(S, model.domain, system, model.formulation)
end

function select_parameters!(S, system::MultiPhaseSystem, model)
    select_default_darcy_parameters!(S, model.domain, system, model.formulation)
end

function select_default_darcy_secondary_variables!(S, domain, system, formulation)
    nph = number_of_phases(system)
    S[:PhaseMassDensities] = ConstantCompressibilityDensities(nph)
    S[:TotalMasses] = TotalMasses()
    if !(isa(system, SinglePhaseSystem) || isa(domain.grid, WellGrid))
        S[:RelativePermeabilities] = BrooksCoreyRelPerm(system)
    end
end

function select_default_darcy_parameters!(prm, domain, system, formulation)
    prm[:PhaseViscosities] = PhaseViscosities() # ConstantVariables(1e-3*ones(nph), Cells()) # 1 cP for all phases by default
    prm[:FluidVolume] = FluidVolume()
    if isa(system, SinglePhaseSystem) || isa(domain.grid, WellGrid)
        prm[:RelativePermeabilities] = RelativePermeabilities()
    end
    if isa(system, SinglePhaseSystem)
        prm[:Saturations] = Saturations()
    end
end

function select_minimum_output_variables!(out, system::MultiPhaseSystem, model)
    push!(out, :TotalMasses)
end

struct BrooksCoreyRelPerm{V, T} <: AbstractRelativePermeabilities
    exponents::V
    residuals::V
    endpoints::V
    residual_total::T
    function BrooksCoreyRelPerm(sys_or_nph::Union{MultiPhaseSystem, Integer}, exponents = 1.0, residuals = 0.0, endpoints = 1.0)
        if isa(sys_or_nph, Integer)
            nph = sys_or_nph
        else
            nph = number_of_phases(sys_or_nph)
        end
        e = expand_to_phases(exponents, nph)
        r = expand_to_phases(residuals, nph)
        epts = expand_to_phases(endpoints, nph)

        total = sum(residuals)
        new{typeof(e), typeof(total)}(e, r, epts, total)
    end
end

function transfer(c::SingleCUDAContext, kr::BrooksCoreyRelPerm)
    e = transfer(c, kr.exponents)
    r = transfer(c, kr.residuals)
    ept = transfer(c, kr.residual_total)

    nph = length(e)
    BrooksCoreyRelPerm(nph, e, r, ept)
end

@jutul_secondary function update_as_secondary!(kr, kr_def::BrooksCoreyRelPerm, model, Saturations)
    n, sr, kwm, sr_tot = kr_def.exponents, kr_def.residuals, kr_def.endpoints, kr_def.residual_total
    @tullio kr[ph, i] = brooks_corey_relperm(Saturations[ph, i], n[ph], sr[ph], kwm[ph], sr_tot)
end

function brooks_corey_relperm(s::T, n::Real, sr::Real, kwm::Real, sr_tot::Real) where T
    den = 1 - sr_tot
    sat = (s - sr) / den
    sat = clamp(sat, zero(T), one(T))
    return kwm*sat^n
end

"""
Interpolated multiphase rel. perm. that is simple (single region, no magic for more than two phases)
"""
struct TabulatedRelPermSimple{V, M, I} <: AbstractRelativePermeabilities
    s::V
    kr::M
    interpolators::I
    function TabulatedRelPermSimple(s::AbstractVector, kr::AbstractVector; kwarg...)
        nph = length(kr)
        n = length(kr[1])
        @assert nph > 0
        T = eltype(kr[1])
        if n <= 50
            V = SVector{n, T}
        else
            V = Vector{T}
        end
        if eltype(s)<:AbstractVector
            # We got a set of different vectors that correspond to rows of kr
            @assert all(map(length, s) .== map(length, kr))
            interpolators = map((ix) -> get_1d_interpolator(V(s[ix]), V(kr[ix]); kwarg...), 1:nph)
        else
            # We got a single vector that is used for all rows
            @assert length(s) == n
            interpolators = map((ix) -> get_1d_interpolator(V(s), V(kr[ix]); kwarg...), 1:nph)
        end
        i_t = Tuple(interpolators)
        new{typeof(s), typeof(kr), typeof(i_t)}(s, kr, i_t)
    end
end

@jutul_secondary function update_as_secondary!(kr, kr_def::TabulatedRelPermSimple, model, Saturations)
    I = kr_def.interpolators
    if false
        @tullio kr[ph, i] = I[ph](Saturations[ph, i])
    else
        tb = minbatch(model.context)
        threaded_interp!(kr, model.context, I, Saturations)
    end
end

function threaded_interp!(F, context, I, x)
    tb = minbatch(context)
    nc = size(F, 2)
    apply(I, x, j, i) = @inbounds I(x[j, i])
    @batch minbatch = tb for i in 1:nc
        @inbounds for j in eachindex(I)
            F[j, i] = apply(I[j], x, j, i)
        end
    end
end

"""
Interpolated multiphase rel. perm. that is simple (single region, no magic for more than two phases)
"""
struct ThreePhaseRelPerm{O, OW, OG, G, R} <: AbstractRelativePermeabilities
    krw::O
    krow::OW
    krog::OG
    krg::G
    swcon::R
end

function ThreePhaseRelPerm(; w, g, ow, og, swcon = 0.0)
    return ThreePhaseRelPerm(w, ow, og, g, swcon)
end

@jutul_secondary function update_as_secondary!(kr, relperm::ThreePhaseRelPerm, model, Saturations)
    s = Saturations
    krw = relperm.krw
    krow = relperm.krow
    krog = relperm.krog
    krg = relperm.krg
    swcon = relperm.swcon

    l = 1
    o = 2
    g = 3
    @inbounds for c in 1:size(kr, 2)
        # Water
        sw = s[l, c]
        kr[l, c] = krw(sw)
        # Gas
        sg = s[g, c]
        kr[g, c] = krg(sg)
        # Oil is special
        so = s[o, c]
        swc = min(swcon, value(sw) - 1e-5)
        d  = (sg + sw - swc)
        ww = (sw - swc)/d
        kro = (1-ww)*krog(so) + ww*krow(so)
        kr[o, c] = kro
    end
end

struct SimpleCapillaryPressure{T} <: JutulVariables
    pc::T
end

degrees_of_freedom_per_entity(model, v::SimpleCapillaryPressure) = number_of_phases(model.system) - 1

@jutul_secondary function update_as_secondary!(Δp, pc::SimpleCapillaryPressure, model, Saturations)
    cap = pc.pc
    npc, nc = size(Δp)
    if npc == 1
        pcow = cap[1]
        @inbounds for c in 1:nc
            sw = Saturations[1, c]
            Δp[1, c] = pcow(sw)
        end
    elseif npc == 2
        pcow, pcog = cap
        if isnothing(pcow)
            @inbounds for c in 1:nc
                sg = Saturations[3, c]
                Δp[1, c] = 0
                Δp[2, c] = pcog(sg)
            end
        elseif isnothing(pcog)
            @inbounds for c in 1:nc
                sw = Saturations[1, c]
                Δp[1, c] = -pcow(sw)
                Δp[2, c] = 0
            end
        else
            @inbounds for c in 1:nc
                sw = Saturations[1, c]
                sg = Saturations[3, c]
                Δp[1, c] = -pcow(sw)
                Δp[2, c] = pcog(sg)
            end
        end
    else
        error("Only implemented for two and three-phase flow.")
    end
end


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

expand_to_phases(v::Real, nph) = SVector{nph}([v for i in 1:nph])
function expand_to_phases(v::AbstractVector, nph)
    @assert length(v) == nph
    return SVector{nph}(v)
end
