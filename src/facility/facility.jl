"""
(Absolute) Minimum well rate for a well that is not disabled.
"""
const MIN_ACTIVE_WELL_RATE = 1e-20
"""
(Absolute) Minimum initial rate for wells when controls are updated.
"""
const MIN_INITIAL_WELL_RATE = 1e-12
"""
Well variables - entities that we have exactly one of per well (and usually relates to the surface connection)
"""

include("types.jl")
include("wells/wells.jl")
include("controls.jl")
include("wellgroups.jl")
include("cross_terms.jl")
include("well_presolve.jl")
include("gradients.jl")

function Jutul.select_minimum_output_variables!(vars, domain::Union{MSWellDomain, SimpleWellDomain}, model::SimulationModel{<:Any, CompositeSystem{:Reservoir, T}, <:Any, <:Any}) where T
    push!(vars, :PhaseMassDensities)
    push!(vars, :Saturations)
    push!(vars, :SurfaceWellConditions)
    return vars
end

function Jutul.select_minimum_output_variables!(vars, domain::WellGroup, model)
    for k in keys(model.primary_variables)
        push!(vars, k)
    end
    push!(vars, :WellGroupConfiguration)
    return vars
end

function setup_injector_control(val::WellTarget, mix; kwarg...)
    return InjectorControl(val, mix; kwarg...)
end

function setup_injector_control(val, type, mix; kwarg...)
    if type isa AbstractString
        type = Symbol(type)
    end
    info = well_target_information(type)
    if ismissing(info.type)
        error("Unknown well target type '$val' - or missing type field in well_target_information.")
    end
    target = info.type(val)
    return InjectorControl(target, mix; kwarg...)
end

function setup_producer_control(val::WellTarget; kwarg...)
    return ProducerControl(val; kwarg...)
end

function setup_producer_control(val::Number, type::Union{String, Symbol}; signed = false, kwarg...)
    if type isa AbstractString
        type = Symbol(type)
    end
    if !signed
        val = abs(val)
    end
    info = well_target_information(type)
    if ismissing(info.type)
        error("Unknown well target type '$val' - or missing type field in well_target_information.")
    end
    target = info.type(val)

    return ProducerControl(target; kwarg...)
end

function setup_disabled_control()
    return DisabledControl()
end
