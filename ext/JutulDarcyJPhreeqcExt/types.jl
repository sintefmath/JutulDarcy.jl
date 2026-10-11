struct PhreeqcTracer <: Tracers.AbstractTracer
    phase_index::Int
    concentration_scale::Float64
end

Tracers.tracer_phase_indices(t::PhreeqcTracer) = (t.phase_index, )
Tracers.tracer_scale(model, t::PhreeqcTracer) = t.concentration_scale

function Jutul.get_dependencies(::PhreeqcTracer, model)
    return JutulDarcy.model_or_domain_is_well(model) ? () : (:PhreeqcAqueousConcentrations, )
end

@inline function Tracers.tracer_phase_concentration(::PhreeqcTracer, model, state,
        index, phase, cell, C = state.TracerConcentrations)
    if JutulDarcy.model_or_domain_is_well(model)
        return C[index, cell]
    else
        return state.PhreeqcAqueousConcentrations[index, cell]
    end
end

function Tracers.tracer_total_mass_flux(tracer::PhreeqcTracer, model, state,
        phase_mass_fluxes, index, upw, C = state.TracerConcentrations, T = eltype(C))
    phase = tracer.phase_index
    q = phase_mass_fluxes[phase]
    concentration = cell -> Tracers.tracer_phase_concentration(tracer, model, state, index, phase, cell, C)
    return Jutul.upwind(upw, concentration, q)*q
end

struct PhreeqcConcentrations <: Jutul.VectorVariables
    signed::Vector{Bool}
    reference::Vector{Float64}
end
Jutul.values_per_entity(model, c::PhreeqcConcentrations) = length(c.signed)

# Excess H and O can be negative. A common nonnegative tracer bound would
# silently destroy alkalinity; impose bounds individually instead.
function Jutul.update_primary_variable!(state, c::PhreeqcConcentrations, key, model, dx, w)
    vals = state[key]
    active = Jutul.active_entities(model.domain, Cells(), for_variables = true)
    function update(local_cell)
        cell = active[local_cell]
        for i in eachindex(c.signed)
            lo = c.reference[i] > 0 ? 1e-12 - c.reference[i] : (c.signed[i] ? nothing : 0.0)
            vals[i, cell] = Jutul.update_value(vals[i, cell], w*dx[i, local_cell],
                nothing, nothing, lo, nothing, nothing)
        end
        return
    end
    Jutul.threaded_loop_minbatch(update, length(active), model.context)
    return vals
end

mutable struct ChemistryWorker
    id::Int
    x::Matrix{Float64}
    trial::Matrix{Float64}
    concentrations::Vector{Float64}
    pressure::Vector{Float64}
    temperature::Vector{Float64}
    reference::Vector{Float64}
    initial_phases::Vector{Int}
    reuse_speciation::Bool
    needs_restore::Bool
    stoichiometry::Matrix{Float64}
    molar_volumes::Vector{Float64}
    selected_output::Vector{Float64}
    saturation_indices::Matrix{Float64}
    y::Matrix{Float64}
    perturbed::Matrix{Float64}
    jacobian::Array{Float64, 3}
    steps::Vector{Float64}
end

mutable struct PhreeqcEquilibrium <: Jutul.VectorVariables
    workers::Vector{ChemistryWorker}
    available::Channel{ChemistryWorker}
    task_key::Symbol
    components::Vector{String}
    molar_masses::Vector{Float64}
    minerals::Vector{String}
    solutions::Vector{UnitRange{Int}}
    signed::Vector{Bool}
    reference::Vector{Float64}
    phase_index::Int
    temperature::Float64
    variable_temperature::Bool
    difference::Symbol
    relative_step::Float64
    batch_size::Int
    cache_input::Matrix{Float64}
    cache_output::Matrix{Float64}
    cache_jacobian::Array{Float64, 3}
    cache_valid::Vector{Bool}
    cache_derivatives::Vector{Bool}
    cache_locks::Vector{ReentrantLock}
    closed::Bool
    active_updates::Int
    lifecycle_lock::ReentrantLock
    lifecycle_condition::Threads.Condition
end

struct PhreeqcImmobileConcentrations <: Jutul.VectorVariables
    ncomponents::Int
end

struct PhreeqcSplitConcentrations <: Jutul.VectorVariables
    rows::UnitRange{Int}
end

struct PhreeqcMineralVolume <: Jutul.ScalarVariable
    row::Int
end

struct PhreeqcReferenceMineralVolume <: Jutul.ScalarVariable
    values::Vector{Float64}
end

struct PhreeqcPoreVolume <: Jutul.ScalarVariable
    base_volume::Symbol
    phase_index::Int
end

Jutul.default_parameter_values(data_domain, model, v::PhreeqcReferenceMineralVolume, key) = copy(v.values)
Jutul.minimum_value(::PhreeqcReferenceMineralVolume) = 0.0
Jutul.minimum_value(::PhreeqcPoreVolume) = eps()

Jutul.values_per_entity(model, c::PhreeqcEquilibrium) = length(c.components) + length(c.minerals) + 1
Jutul.values_per_entity(model, c::PhreeqcImmobileConcentrations) = c.ncomponents
Jutul.values_per_entity(model, c::PhreeqcSplitConcentrations) = length(c.rows)

# Native handles must never be duplicated by copying Julia structs. Model and
# storage copies share this synchronized runtime; cached values contain no AD.
Base.deepcopy_internal(c::PhreeqcEquilibrium, dict::IdDict) = c

function Jutul.subvariable(c::PhreeqcEquilibrium, map::Jutul.FiniteVolumeGlobalMap)
    error("PHREEQC chemistry requires a separate setup for each subdomain; NLDD/MPI partitioning is not supported.")
end
