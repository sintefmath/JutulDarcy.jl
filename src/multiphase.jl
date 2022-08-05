export MultiPhaseSystem, ImmiscibleSystem, SinglePhaseSystem
export AqueousPhase, LiquidPhase, VaporPhase
export number_of_phases, get_short_name, get_name, subscript
export update_linearized_system!
export SourceTerm
export setup_storage, update_equations!
export Pressure, Saturations, TotalMasses, TotalMass

# Abstract multiphase system

get_phases(sys::MultiPhaseSystem) = sys.phases
number_of_phases(sys::MultiPhaseSystem) = length(get_phases(sys))

@enum FlowSourceType begin
    MassSource
    StandardVolumeSource
    VolumeSource
end

struct SourceTerm{I, F, T} <: JutulForce
    cell::I
    value::F
    fractional_flow::T
    type::FlowSourceType
    function SourceTerm(cell, value; fractional_flow = [1.0], type = MassSource)
        @assert sum(fractional_flow) == 1.0 "Fractional flow for source term in cell $cell must sum to 1."
        f = Tuple(fractional_flow)
        return new{typeof(cell), typeof(value), typeof(f)}(cell, value, f, type)
    end
end

function cell(s::SourceTerm{I, T}) where {I, T} 
    return s.cell::I
end

function setup_forces(model::SimulationModel{G, S}; sources = nothing) where {G<:Any, S<:MultiPhaseSystem}
    return (sources = sources,)
end

function Jutul.subforce(s::AbstractVector{S}, model) where S<:SourceTerm
    # Just to be safe
    s = deepcopy(s)
    m = global_map(model.domain)

    n = length(s)
    keep = repeat([false], n)
    for (i, src) in enumerate(s)
        # Cell must be in local domain, and not on boundary
        if !global_cell_inside_domain(src.cell, m)
            continue
        end
        c_l = local_cell(src.cell, m)
        c_i = interior_cell(c_l, m)
        inner = !isnothing(c_i)
        if !inner
            continue
        end
        keep[i] = true
        s[i] = SourceTerm(c_i, src.value, fractional_flow = src.fractional_flow)
    end
    return s[keep]
end


number_of_components(sys::ImmiscibleSystem) = number_of_phases(sys)

# function ImmiscibleSystem(phases)
#    @assert length(phases) > 1 "System should have at least two phases. For single-phase, use SinglePhaseSystem instead."
# end

# Single-phase

phase_names(system) = get_name.(get_phases(system))

get_phases(sys::SinglePhaseSystem) = (sys.phase, )
number_of_phases(::SinglePhaseSystem) = 1

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

"""
Pressure
"""
struct Pressure <: ScalarVariable
    max_abs
    max_rel
    minimum_pressure
    maximum_pressure
    scale
    function Pressure(; max_abs = nothing, max_rel = nothing, scale = 1e8, maximum = Inf, minimum = -Inf)
        new(max_abs, max_rel, minimum, maximum, scale)
    end
end

Jutul.variable_scale(p::Pressure) = p.scale
absolute_increment_limit(p::Pressure) = p.max_abs
relative_increment_limit(p::Pressure) = p.max_rel
maximum_value(p::Pressure) = p.maximum_pressure
minimum_value(p::Pressure) = p.minimum_pressure

# Saturations as primary variable
struct Saturations <: FractionVariables
    ds_max
    Saturations(;ds_max = 0.2) = new(ds_max)
end

values_per_entity(model, v::Saturations) = number_of_phases(model.system)
absolute_increment_limit(s::Saturations) = s.ds_max

function initialize_primary_variable_ad!(state, model, pvar::Saturations, state_symbol, npartials; kwarg...)
    nph = values_per_entity(model, pvar)
    v = state[state_symbol]
    state[state_symbol] = Jutul.unit_sum_init(v, model, npartials, nph; kwarg...)
    return state
end

# Total component masses
struct TotalMasses <: GroupedVariables end

function degrees_of_freedom_per_entity(model::SimulationModel{G, S}, v::TotalMasses) where {G<:Any, S<:MultiPhaseSystem}
    number_of_phases(model.system)
end

@inline function minimum_value(::TotalMasses) 0 end

struct TotalMass <: ScalarVariable end
@inline function minimum_value(::TotalMass) 0 end


# Selection of variables
function select_primary_variables_system!(S, domain, system::SinglePhaseSystem, formulation)
    S[:Pressure] = Pressure()
end

function select_primary_variables_system!(S, domain, system::ImmiscibleSystem, formulation)
    S[:Pressure] = Pressure()
    S[:Saturations] = Saturations()
end

function select_equations_system!(eqs, domain, system::MultiPhaseSystem, formulation)
    eqs[:mass_conservation] = ConservationLaw(domain.discretizations.mass_flow)
end

number_of_equations_per_entity(system::MultiPhaseSystem, e::ConservationLaw) = number_of_components(system)
number_of_equations_per_entity(system::SinglePhaseSystem, e::ConservationLaw) = 1

export fluid_volume, pore_volume
pore_volume(model::MultiModel) = pore_volume(reservoir_model(model))
pore_volume(model::SimulationModel) = fluid_volume(model.domain.grid)
pore_volume(grid) = fluid_volume(grid)

fluid_volume(domain::DiscretizedDomain) = fluid_volume(domain.grid)
fluid_volume(grid::MinimalTPFAGrid) = grid.pore_volumes
fluid_volume(grid) = 1.0


function apply_forces_to_equation!(storage, model::SimulationModel{D, S}, eq::ConservationLaw, force::V, time) where {V <: AbstractVector{SourceTerm{I, F, T}}, D, S<:MultiPhaseSystem} where {I, F, T}
    acc = get_diagonal_entries(eq)
    state = storage.state
    if haskey(state, :RelativePermeabilities)
        kr = state.RelativePermeabilities
    else
        kr = 1.0
    end
    mu = state.PhaseViscosities
    rhoS = get_reference_densities(model, storage)
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

function insert_phase_sources!(acc::CuArray, model, kr, mu, rhoS, sources)
    sources::CuArray
    ix = map(cell, sources)
    if !isa(rhoS, CuArray)
        @warn "SurfaceDensities is not a CuArray, will convert whenever needed. Improve performance by converting once."  maxlog=1
        rhoS = CuArray(rhoS)
    end
    rhoS::CuArray
    @tullio acc[ph, ix[i]] = acc[ph, ix[i]] - phase_source(sources[i].cell, sources[i], rhoS[ph], kr, mu, ph)
end

function convergence_criterion(model::SimulationModel{D, S}, storage, eq::ConservationLaw, eq_s, r; dt = 1) where {D, S<:MultiPhaseSystem}
    M = global_map(model.domain)
    v = x -> Jutul.active_view(x, M, for_variables = false)
    Φ = v(storage.state.FluidVolume)
    ρ = v(storage.state.PhaseMassDensities)

    @timeit "cnv" @tullio max e[j] := abs(r[j, i]) * dt / (value(ρ[j, i])*value(Φ[i]))
    @timeit "mb" begin
        scale = dt./(sum(value, Φ).*mean(value, ρ, dims = 2))
        mb = scale.*sum(abs, r, dims = 2)
    end

    names = phase_names(model.system)
    R = Dict("CNV" => (errors = e, names = names),
             "MB"  => (errors = mb, names = names))
    return R
end

function get_reference_densities(model, storage)
    prm = storage.parameters
    return typed_reference_density(prm.reference_densities, model)
end

typed_reference_density(rhoS::Tuple, model) = rhoS
function typed_reference_density(ρ, model)
    N = length(ρ)
    T = Jutul.float_type(model.context)
    return tuple(ρ...)::NTuple{N, T}
end

function setup_parameters_system!(d, model, sys::MultiPhaseSystem)
    nph = number_of_phases(sys)
    d[:reference_densities] = transfer(model.context, ones(nph))
end

function cpr_weights_no_partials!(w, model::SimulationModel{R, S}, state, r, n, bz, scaling) where {R, S<:ImmiscibleSystem}
    ρ = state.PhaseMassDensities
    tb = minbatch(model.context)
    nph, nc = size(ρ)
    @batch minbatch = tb for i in 1:nc
        for ph in 1:nph
            @inbounds w[ph, i] = 1/value(ρ[ph, i])
        end
    end
end

function capillary_pressure(model, s)
    pck = :CapillaryPressures
    if haskey(s, pck)
        pc = s[pck]
        ref_index = min(2, number_of_phases(model.system))
    else
        pc = nothing
        ref_index = 1
    end
    return(pc, ref_index)
end
