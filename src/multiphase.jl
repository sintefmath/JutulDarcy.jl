export MultiPhaseSystem, ImmiscibleSystem, SinglePhaseSystem
export AqueousPhase, LiquidPhase, VaporPhase
export number_of_phases, get_short_name, get_name, subscript
export update_linearized_system!
export SourceTerm
export setup_storage, update_equations!
export Pressure, Saturations, TotalMasses, TotalMass

# Abstract multiphase system

get_phases(sys::MultiPhaseSystem) = sys.phases
@inline number_of_phases(sys::MultiPhaseSystem) = length(get_phases(sys))
@inline reference_densities(sys::MultiPhaseSystem) = sys.rho_ref
@inline reference_densities(sys::CompositeSystem) = reference_densities(flow_system(sys))

flow_system(sys::MultiPhaseSystem) = sys
flow_system(sys::CompositeSystem) = sys.systems.flow





number_of_components(sys::ImmiscibleSystem) = number_of_phases(sys)

# function ImmiscibleSystem(phases)
#    @assert length(phases) > 1 "System should have at least two phases. For single-phase, use SinglePhaseSystem instead."
# end

# Single-phase

phase_names(system) = get_name.(get_phases(system))

get_phases(sys::SinglePhaseSystem) = (sys.phase, )

phase_indices(sys::SinglePhaseSystem) = 1
phase_indices(sys::ImmiscibleSystem) = tuple(eachindex(sys.phases)...)


number_of_phases(::SinglePhaseSystem) = 1
number_of_phases(sys::CompositeSystem) = number_of_phases(sys.systems.flow)

## Phases
# Abstract phase
abstract type AbstractPhase end

function get_short_name(phase::AbstractPhase)
    return get_name(phase)[]
end

function subscript(prefix::String, phase::AbstractPhase)
    return string(prefix, "_", get_short_name(phase))
end
# Aqueous phase
struct AqueousPhase <: AbstractPhase end
get_name(::AqueousPhase) = "Aqueous"

# Liquid phase
struct LiquidPhase <: AbstractPhase end
get_name(::LiquidPhase) = "Liquid"

# Vapor phases
struct VaporPhase <: AbstractPhase end
get_name(::VaporPhase) = "Vapor"

## Main implementation
# Primary variable logic

const DEFAULT_MINIMUM_PRESSURE = 101325.0

struct Pressure <: ScalarVariable
    max_abs::Union{Float64, Nothing}
    max_rel::Union{Float64, Nothing}
    minimum_pressure::Float64
    maximum_pressure::Float64
    scale::Float64
end

"""
    Pressure(; max_abs = nothing, max_rel = 0.2, scale = 1e8, maximum = Inf, minimum = DEFAULT_MINIMUM_PRESSURE)

Pressure variable definition. `max_abs`/`max_rel` maximum allowable
absolute/relative change over a Newton iteration, `scale` is a "typical" value
used to regularize the linear system, `maximum` the largest possible value and
`minimum` the smallest.
"""
function Pressure(; max_abs = nothing, max_rel = 0.2, scale = 1e8, maximum = Inf, minimum = DEFAULT_MINIMUM_PRESSURE)
    Pressure(max_abs, max_rel, minimum, maximum, scale)
end

Jutul.variable_scale(p::Pressure) = p.scale
absolute_increment_limit(p::Pressure) = p.max_abs
relative_increment_limit(p::Pressure) = p.max_rel
maximum_value(p::Pressure) = p.maximum_pressure
minimum_value(p::Pressure) = p.minimum_pressure

struct Saturations <: FractionVariables
    ds_max::Float64
end

"""
    Saturations(;ds_max = 0.2)

Saturations as primary variable. `ds_max` controls maximum allowable saturation
change between two Newton iterations.
"""
function Saturations(;ds_max = 0.2)
    Saturations(ds_max)
end

default_value(model::SimulationModel{<:Any, <:SinglePhaseSystem, <:Any, <:Any}, ::Saturations) = 1.0

values_per_entity(model, v::Saturations) = number_of_phases(model.system)
absolute_increment_limit(s::Saturations) = s.ds_max

function initialize_primary_variable_ad!(state, model, pvar::Saturations, state_symbol, npartials; kwarg...)
    nph = values_per_entity(model, pvar)
    v = state[state_symbol]
    state[state_symbol] = Jutul.unit_sum_init(v, model, npartials, nph; kwarg...)
    return state
end

# Total component masses
struct TotalMasses <: VectorVariables end

function degrees_of_freedom_per_entity(model::SimulationModel{G, S}, v::TotalMasses) where {G<:Any, S<:MultiPhaseSystem}
    number_of_phases(model.system)
end

struct PhaseMassMobilities <: PhaseVariables end
struct PhaseMobilities <: PhaseVariables end

@inline function minimum_value(::TotalMasses) 0 end

struct TotalMass <: ScalarVariable end
@inline function minimum_value(::TotalMass) 0 end

struct Transmissibilities <: ScalarVariable end
Jutul.variable_scale(::Transmissibilities) = 1e-10
Jutul.minimum_value(::Transmissibilities) = 0.0

Jutul.associated_entity(::Transmissibilities) = Faces()
function Jutul.default_values(model, ::Transmissibilities)
    return physical_representation(model.domain).trans
end

function Jutul.default_parameter_values(data_domain, model, param::Transmissibilities, symb)
    if haskey(data_domain, :transmissibilities, Faces())
        # This takes precedence
        T = data_domain[:transmissibilities]
    elseif haskey(data_domain, :permeability, Cells())
        U = data_domain[:permeability]
        g = physical_representation(data_domain)
        T = compute_face_trans(g, U)
        if any(x -> x < 0, T)
            c = count(x -> x < 0, T)
            @warn "$c negative transmissibilities detected."
        end
    else
        error(":permeability or :transmissibilities symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

struct TwoPointGravityDifference <: ScalarVariable end

Jutul.associated_entity(::TwoPointGravityDifference) = Faces()
function Jutul.default_values(model, ::TwoPointGravityDifference)
    return physical_representation(model.domain).gdz
end

function Jutul.default_parameter_values(data_domain, model, param::TwoPointGravityDifference, symb)
    has_gdz = haskey(data_domain, :gdz, Faces())
    has_cell_centroids = haskey(data_domain, :cell_centroids, Cells())
    has_face_neighbors = haskey(data_domain, :neighbors, Faces())
    if has_gdz
        # This takes precedence
        gdz = data_domain[:gdz]
    elseif has_cell_centroids && has_face_neighbors
        N = data_domain[:neighbors]
        cc = data_domain[:cell_centroids]
        if size(cc, 1) == 3
            z = vec(cc[3, :])
            gdz = compute_face_gdz(N, z)
        else
            nf = size(N, 2)
            gdz = zeros(nf)
        end
    else
        error(":gdz or :neighbors + :cell_centroids symbols must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return gdz
end

# Selection of variables
function select_primary_variables!(S, ::SinglePhaseSystem, model::SimulationModel)
    S[:Pressure] = Pressure()
end

function select_primary_variables!(S, ::ImmiscibleSystem, model::SimulationModel)
    S[:Pressure] = Pressure()
    S[:Saturations] = Saturations()
end

function select_equations!(eqs, sys::MultiPhaseSystem, model::SimulationModel)
    fdisc = model.domain.discretizations.mass_flow
    nc = number_of_components(sys)
    eqs[:mass_conservation] = ConservationLaw(fdisc, :TotalMasses, nc)
end

function select_parameters!(prm, disc::D, model::DarcyFlowModel) where D<:Union{TwoPointPotentialFlowHardCoded, Jutul.PotentialFlow}
    prm[:Transmissibilities] = Transmissibilities()
    prm[:TwoPointGravityDifference] = TwoPointGravityDifference()
end

number_of_equations_per_entity(system::MultiPhaseSystem, e::ConservationLaw) = number_of_components(system)
number_of_equations_per_entity(system::SinglePhaseSystem, e::ConservationLaw) = 1

export fluid_volume, pore_volume
function pore_volume(data_domain::DataDomain)
    if haskey(data_domain, :pore_volume, Cells())
        vol = data_domain[:pore_volume]
    elseif haskey(data_domain, :volumes, Cells())
        vol = data_domain[:volumes]
        if haskey(data_domain, :porosity, Cells())
            vol = vol.*data_domain[:porosity]
        end
        if haskey(data_domain, :net_to_gross, Cells())
            vol = vol.*data_domain[:net_to_gross]
        end
    else
        error("Neither pair :volumes and :porosity or :pore_volume found in domain.")
    end

    return vol
end

pore_volume(model::MultiModel, parameters) = pore_volume(reservoir_model(model), parameters[:Reservoir])
pore_volume(model::SimulationModel, parameters) = fluid_volume(model, parameters)
fluid_volume(model, parameters) = parameters[:FluidVolume]
domain_fluid_volume(g) = missing

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw, eq_s, force::V, time) where {V <: AbstractVector{SourceTerm{I, F, T}}, D, S<:MultiPhaseSystem} where {I, F, T}
    state = storage.state
    if haskey(state, :RelativePermeabilities)
        kr = state.RelativePermeabilities
    else
        kr = 1.0
    end
    mu = state.PhaseViscosities
    rhoS = reference_densities(model.system)
    insert_phase_sources!(acc, model, kr, mu, rhoS, force)
end

function local_mobility(kr::Real, mu, ph, c)
    return kr./mu[ph, c]
end

function local_mobility(kr, mu, ph, c)
    return kr[ph, c]./mu[ph, c]
end

function phase_source(c, src, rhoS, kr, mu, ph)
    v = src.value
    # c = src.cell
    if v > 0
        q = in_phase_source(src, v, c, kr, mu, ph)
    else
        q = out_phase_source(src, v, c, kr, mu, ph)
    end
    return rhoS*q
end


function in_phase_source(src, v, c, kr, mu, ph)
    f = src.fractional_flow[ph]
    return v*f
end

function out_phase_source(src, v, c, kr, mu, ph)
    mobT = 0
    mob = 0
    for i = 1:size(mu, 1)
        mi = local_mobility(kr, mu, i, c)
        mobT += mi
        if ph == i
            mob = mi
        end
    end
    f = mob/mobT
    return v*f
end

function insert_phase_sources!(acc, model, kr, mu, rhoS, sources)
    nph = size(acc, 1)
    M = global_map(model.domain)
    for src in sources
        c = Jutul.full_cell(src.cell, M)
        for ph = 1:nph
            q_ph = phase_source(c, src, rhoS[ph], kr, mu, ph)
            @inbounds acc[ph, src.cell] -= q_ph
        end
    end
end

# function insert_phase_sources!(acc::CuArray, model, kr, mu, rhoS, sources)
#     sources::CuArray
#     ix = map(cell, sources)
#     if !isa(rhoS, CuArray)
#         @warn "SurfaceDensities is not a CuArray, will convert whenever needed. Improve performance by converting once."  maxlog=1
#         rhoS = CuArray(rhoS)
#     end
#     rhoS::CuArray
#     @tullio acc[ph, ix[i]] = acc[ph, ix[i]] - phase_source(sources[i].cell, sources[i], rhoS[ph], kr, mu, ph)
# end

function convergence_criterion(model::SimulationModel{D, S}, storage, eq::ConservationLaw, eq_s, r; dt = 1, update_report = missing) where {D, S<:MultiPhaseSystem}
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = v(storage.state.FluidVolume)
    ρ = v(storage.state.PhaseMassDensities)

    nph = number_of_phases(model.system)
    cnv, mb = cnv_mb_errors(r, Φ, ρ, dt, Val(nph))

    names = phase_names(model.system)
    R = (CNV = (errors = cnv, names = names),
         MB = (errors = mb, names = names))
    return R
end

function cnv_mb_errors(r, Φ, ρ, dt, ::Val{N}) where N
    nc = length(Φ)
    mb = @MVector zeros(N)
    cnv = @MVector zeros(N)
    avg_density = @MVector zeros(N)

    pv_t = 0.0
    @inbounds for c in 1:nc
        pv_c = Φ[c]
        pv_t += pv_c
        @inbounds for ph = 1:N
            r_ph = r[ph, c]
            ρ_ph = ρ[ph, c]
            # MB
            mb[ph] += r_ph
            avg_density[ph] += ρ_ph
            # CNV
            cnv[ph] = max(cnv[ph], dt*abs(r_ph)/(ρ_ph*pv_c))
        end
    end
    @inbounds for ph = 1:N
        ρ_avg = avg_density[ph]/nc
        mb[ph] = (dt/pv_t)*abs(mb[ph])/ρ_avg
    end
    return (Tuple(cnv), Tuple(mb))
end

function cpr_weights_no_partials!(w, model::SimulationModel{R, S}, state, r, n, bz, scaling) where {R, S<:ImmiscibleSystem}
    ρ = state.PhaseMassDensities
    nc = size(w, 2)
    tb = minbatch(model.context, nc)
    M = global_map(model.domain)
    density = Jutul.active_view(ρ, M, for_variables = false)
    @batch minbatch = tb for i in axes(w, 2)
        for ph in axes(w, 1)
            @inbounds w[ph, i] = 1/value(density[ph, i])
        end
    end
end

function capillary_pressure(model, s)
    pck = :CapillaryPressure
    if haskey(s, pck)
        pc = s[pck]
        ref_index = min(2, number_of_phases(model.system))
    else
        pc = nothing
        ref_index = 1
    end
    return(pc, ref_index)
end
