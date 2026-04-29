"""
    SubdomainPropertyEvaluator

Holds a one-cell submodel and its pre-allocated parameters, used to evaluate
secondary variables (phase densities, enthalpies, etc.) at boundary / injection
conditions without running the full reservoir EOS on every Newton iteration.
"""
struct SubdomainPropertyEvaluator{M, P}
    model::M
    parameters::P
end

"""
    setup_bc_property_evaluator(reservoir_model, cell)

Build a [`SubdomainPropertyEvaluator`](@ref) for a single reservoir `cell`.
"""
function setup_bc_property_evaluator(rmodel::SimulationModel, cell::Int)
    sub = submodel(rmodel, [cell])
    prm = setup_parameters(sub)
    return SubdomainPropertyEvaluator(sub, prm)
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _bc_state_dict(sub_model::SimulationModel, p, T, saturations_or_mix)
    pvars = sub_model.primary_variables
    state = Dict{Symbol, Any}()
    state[:Pressure] = [p]
    if haskey(pvars, :Temperature)
        state[:Temperature] = [T]
    end
    if haskey(pvars, :Saturations)
        nph = length(saturations_or_mix)
        s = reshape(collect(Float64, saturations_or_mix), nph, 1)
        state[:Saturations] = s
    end
    if haskey(pvars, :ImmiscibleSaturation)
        state[:ImmiscibleSaturation] = [Float64(first(saturations_or_mix))]
    end
    if haskey(pvars, :OverallMoleFractions)
        ncomp = length(saturations_or_mix)
        z = reshape(collect(Float64, saturations_or_mix), ncomp, 1)
        state[:OverallMoleFractions] = z
    end
    return state
end

"""
    evaluate_bc_state(ev, p, T, saturations_or_mix)

Run the EOS for a single cell, returning a `Dict{Symbol, Any}` that contains
all secondary variables (e.g. `:PhaseMassDensities`, `:FluidEnthalpy`, …)
evaluated at the given pressure `p`, temperature `T`, and
saturations / overall mole fractions `saturations_or_mix`.
"""
function evaluate_bc_state(ev::SubdomainPropertyEvaluator, p, T, saturations_or_mix)
    state = _bc_state_dict(ev.model, p, T, saturations_or_mix)
    # For compositional models Temperature is a parameter, not a primary variable.
    prm = ev.parameters
    if haskey(prm, :Temperature)
        prm = copy(prm)          # shallow copy — safe to overwrite scalar entry
        prm[:Temperature] = [Float64(T)]
    end
    return evaluate_all_secondary_variables(ev.model, state, prm)
end

# ---------------------------------------------------------------------------
# Single-force overloads
# ---------------------------------------------------------------------------

"""
    with_property_evaluator(reservoir_model, bc::FlowBoundaryCondition)

Return a new [`FlowBoundaryCondition`](@ref) with the `state_bc` field
populated by evaluating the reservoir EOS at `(bc.pressure, bc.temperature)`
using a one-cell submodel built at `bc.cell`.
"""
function with_property_evaluator(rmodel::SimulationModel, bc::FlowBoundaryCondition)
    nph = number_of_phases(rmodel.system)
    f_inj = bc.fractional_flow
    saturations = if isnothing(f_inj)
        fill(1.0 / nph, nph)
    else
        collect(f_inj)
    end
    ev = setup_bc_property_evaluator(rmodel, bc.cell)
    state_bc = evaluate_bc_state(ev, bc.pressure, bc.temperature, saturations)
    return FlowBoundaryCondition(
        bc.cell, bc.pressure, bc.temperature,
        bc.trans_flow, bc.trans_thermal,
        bc.fractional_flow, bc.density,
        state_bc
    )
end

"""
    with_property_evaluator(reservoir_model, ctrl::InjectorControl, perf_cell)

Return a new [`InjectorControl`](@ref) with the `state_well` field populated
by evaluating the reservoir EOS at `(p_ref, ctrl.temperature)` using a one-cell
submodel built at `perf_cell`.  `p_ref` defaults to `DEFAULT_MINIMUM_PRESSURE`.
"""
function with_property_evaluator(rmodel::SimulationModel, ctrl::InjectorControl, perf_cell::Int; p_ref = DEFAULT_MINIMUM_PRESSURE)
    ev = setup_bc_property_evaluator(rmodel, perf_cell)
    state_well = evaluate_bc_state(ev, p_ref, ctrl.temperature, ctrl.injection_mixture)
    return InjectorControl(
        ctrl.target, ctrl.injection_mixture;
        density = ctrl.mixture_density,
        phases = ctrl.phases,
        temperature = ctrl.temperature,
        enthalpy = ctrl.enthalpy,
        tracers = ctrl.tracers,
        factor = ctrl.factor,
        state_well = state_well,
        check = false
    )
end

# ---------------------------------------------------------------------------
# Recursive multi-force overloads
# ---------------------------------------------------------------------------

with_property_evaluators(model, ::Nothing) = nothing

"""
    with_property_evaluators(reservoir_model, bcs::AbstractVector{<:FlowBoundaryCondition})

Apply [`with_property_evaluator`](@ref) to every boundary condition in `bcs`.
"""
function with_property_evaluators(rmodel::SimulationModel, bcs::AbstractVector{<:FlowBoundaryCondition})
    return map(bc -> with_property_evaluator(rmodel, bc), bcs)
end

# Reservoir sub-forces NamedTuple  (bc = ..., sources = ...)
function _with_property_evaluators_reservoir(rmodel, forces::NamedTuple)
    bc = get(forces, :bc, nothing)
    new_bc = with_property_evaluators(rmodel, bc)
    return merge(forces, (bc = new_bc,))
end

# Facility sub-forces NamedTuple (control = ..., limits = ...)
function _with_property_evaluators_facility(multimodel::MultiModel, rmodel, forces::NamedTuple)
    control = get(forces, :control, nothing)
    isnothing(control) && return forces
    new_control = copy(control)
    for (wname, ctrl) in pairs(control)
        ctrl isa InjectorControl || continue
        haskey(multimodel.models, wname) || continue
        well_domain = multimodel.models[wname].domain
        perf_cell = physical_representation(well_domain).perforations.reservoir[1]
        new_control[wname] = with_property_evaluator(rmodel, ctrl, perf_cell)
    end
    return merge(forces, (control = new_control,))
end

"""
    with_property_evaluators(model::MultiModel, forces::NamedTuple)

Process a full forces `NamedTuple` as returned by `setup_reservoir_forces`,
returning a new `NamedTuple` where every [`FlowBoundaryCondition`](@ref) has
`state_bc` populated and every [`InjectorControl`](@ref) has `state_well`
populated.

Call this once before running a simulation to avoid repeated EOS evaluations
during Newton iterations.
"""
function with_property_evaluators(model::MultiModel, forces::NamedTuple)
    rmodel = reservoir_model(model)
    result = forces
    if haskey(forces, :Reservoir) && !isnothing(forces.Reservoir)
        result = merge(result, (Reservoir = _with_property_evaluators_reservoir(rmodel, forces.Reservoir),))
    end
    if haskey(forces, :Facility) && !isnothing(forces.Facility)
        result = merge(result, (Facility = _with_property_evaluators_facility(model, rmodel, forces.Facility),))
    end
    return result
end

function with_property_evaluators(model::MultiModel, forces::AbstractDict)
    rmodel = reservoir_model(model)
    result = copy(forces)
    if haskey(forces, :Reservoir) && !isnothing(forces[:Reservoir])
        result[:Reservoir] = _with_property_evaluators_reservoir(rmodel, forces[:Reservoir])
    end
    if haskey(forces, :Facility) && !isnothing(forces[:Facility])
        result[:Facility] = _with_property_evaluators_facility(model, rmodel, forces[:Facility])
    end
    return result
end

"""
    with_property_evaluators(model, forces::AbstractVector)

Apply [`with_property_evaluators`](@ref) element-wise to a vector of forces
(e.g. one entry per time step).
"""
function with_property_evaluators(model, forces::AbstractVector)
    return map(f -> with_property_evaluators(model, f), forces)
end
