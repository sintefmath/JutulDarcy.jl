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

"""
    ctrl = setup_injector_control(val, type, mixture)
    ctrl = setup_injector_control(100.0*si_unit("meter")^3/si_unit("day"), "rate", [0.3, 0.7])

Set up an injector control with a given value, type and mixture. This is a
convenience constructor for [`InjectorControl`](@ref). The `type` can be one of
the following:
- "rate" - total injection rate (m³/s)
- "bhp" - bottom hole pressure (Pa)

# Keyword arguments

- `density`: Specify the surface density of the injected fluid (kg/m³). This
  defaults to 1.0, which may be very low for many mixtures (water being around
  1000 kg/m³).
- `enthalpy`: Specify the injected specific enthalpy.

For more details, see [`InjectorControl`](@ref).

# Notes
The mixture is a vector of fractions for each phase, e.g. [0.3, 0.7] for 30%
water and 70% oil at surface conditions. The mixture must sum to one.
"""
function setup_injector_control(val, type, mix; kwarg...)
    if type isa AbstractString
        type = lowercase(type)
        if type == "wrat" || type == "grat" || type == "orat"
            error("To specify injection mixture, set type to \"rate\" and mixture values based on the phase ordering (if you have water, oil and gas, then water injection is [1.0, 0.0, 0.0]).")
        end
        if type == "wbhp"
            type = :bhp
        end
        type = Symbol(lowercase(type))
    end
    info = well_target_information(type)
    if ismissing(info) || ismissing(info.type)
        error("Unknown well target type '$type' - or missing type field in well_target_information.")
    end
    target = info.type(val)
    return InjectorControl(target, mix; kwarg...)
end

function setup_producer_control(val::WellTarget; kwarg...)
    return ProducerControl(val; kwarg...)
end

"""
    ctrl = setup_producer_control(val, type)
    ctrl = setup_producer_control(100.0*si_unit("bar"), "bhp")
    ctrl = setup_producer_control(100.0*si_unit("meter")^3/si_unit("day"), "orat")

Set up a producer control with a given value and type. This is a convenience
constructor for [`ProducerControl`](@ref). The `type` can be one of the
following:
- "rate" - total production rate (m³/s)
- "bhp" - bottom hole pressure (Pa)
_ "orat" - oil production rate (m³/s)
- "wrat" - water production rate (m³/s)
- "grat" - gas production rate (m³/s)
- "lrat" - liquid production rate (m³/s)
- "wcut" - water cut (fraction, can lead to convergence issues)
- "gor" - gas-oil ratio (m³/m³, can lead to convergence issues)
- "wgr" - water-gas ratio (m³/m³, can lead to convergence issues)
- "glr" - gas-liquid ratio (m³/m³, can lead to convergence issues)
"""
function setup_producer_control(val::Number, type; signed = false, kwarg...)
    if type isa AbstractString
        s = lowercase(type)
        if s == "wbhp"
            s = :bhp
        end
        type = Symbol(s)
    end
    type::Symbol
    info = well_target_information(type)
    if ismissing(info) || ismissing(info.type)
        error("Unknown well target type '$type' - or missing type field in well_target_information.")
    end
    if info.is_rate && !signed
        val = -abs(val)
    end
    target = info.type(val)

    return ProducerControl(target; kwarg...)
end

function setup_disabled_control()
    return DisabledControl()
end
