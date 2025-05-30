function Jutul.active_entities(wg::WellGroup, ::Jutul.TrivialGlobalMap, e; for_variables = false)
    return 1:count_entities(wg, e)
end

function Jutul.count_entities(wg::WellGroup, ::Wells)
    return length(wg.well_symbols)
end

function get_well_position(d, symbol)
    return findfirst(isequal(symbol), d.well_symbols)
end

function Jutul.associated_entity(::TotalSurfaceMassRate) Wells() end
function Jutul.associated_entity(::SurfaceTemperature) Wells() end

function Jutul.update_primary_variable!(state, massrate::TotalSurfaceMassRate, state_symbol, model, dx, w)
    v = state[state_symbol]
    symbols = model.domain.well_symbols
    cfg = state.WellGroupConfiguration
    # Injectors can only have strictly positive injection rates,
    # producers can only have strictly negative and disabled controls give zero rate.
    abs_max = Jutul.absolute_increment_limit(massrate)
    rel_max = Jutul.relative_increment_limit(massrate)
    function do_update!(wcfg, s, v, dx, ctrl)
        return Jutul.update_value(v, dx)
    end
    function do_update!(wcfg, s, v, dx, ctrl::InjectorControl)
        limit_rate = MIN_INITIAL_WELL_RATE
        if v <= limit_rate && v + dx < limit_rate && model.domain.can_shut_injectors
            @debug "$s Approaching zero rate, disabling injector for current step"
            wcfg.operating_controls[s] = DisabledControl()
            next = Jutul.update_value(v, -value(v))
        else
            next = Jutul.update_value(v, dx, abs_max, rel_max, limit_rate, nothing)
        end
        return next
    end
    function do_update!(wcfg, s, v, dx, ctrl::ProducerControl)
        # A significant negative rate is the valid producer control
        limit_rate = -MIN_INITIAL_WELL_RATE
        if v >= limit_rate && v + dx > limit_rate && model.domain.can_shut_producers
            @debug "$s Approaching zero rate, disabling producer for current step"
            wcfg.operating_controls[s] = DisabledControl()
            next = Jutul.update_value(v, -value(v))
        else
            next = Jutul.update_value(v, dx, abs_max, rel_max, nothing, limit_rate)
        end
        return next
    end
    function do_update!(wcfg, s, v, dx, ctrl::DisabledControl)
        # Set value to zero since we know it is correct.
        return Jutul.update_value(v, -value(v))
    end
    @inbounds for i in eachindex(v)
        s = symbols[i]
        v[i] = do_update!(cfg, s, v[i], w*dx[i], operating_control(cfg, s))
    end
end


rate_weighted(t) = true
rate_weighted(::BottomHolePressureTarget) = false
rate_weighted(::DisabledTarget) = false

target_scaling(::Any) = 1.0
target_scaling(::BottomHolePressureTarget) = 1e5
target_scaling(::SurfaceVolumeTarget) = 1.0/(3600*24.0) # weight by day

Jutul.associated_entity(::ControlEquationWell) = Wells()
Jutul.local_discretization(::ControlEquationWell, i) = nothing
function Jutul.update_equation_in_entity!(v, i, state, state0, eq::ControlEquationWell, model, dt, ldisc = local_discretization(eq, i))
    # Set to zero, do actual control via cross terms
    v[] = 0*state.TotalSurfaceMassRate[i]
end

Jutul.associated_entity(::SurfaceTemperatureEquation) = Wells()
Jutul.local_discretization(::SurfaceTemperatureEquation, i) = nothing
function Jutul.update_equation_in_entity!(v, i, state, state0, eq::SurfaceTemperatureEquation, model, dt, ldisc = local_discretization(eq, i))
    # Set equal to surface temperature. corresponding well temperatures will be
    # subtracted using corss terms
    v[] = state.SurfaceTemperature[i]
end

# Selection of primary variables
function select_primary_variables!(S, domain::WellGroup, model)
    S[:TotalSurfaceMassRate] = TotalSurfaceMassRate()
end

function select_primary_variables!(S, system::PredictionMode, model)
    nothing
end

function select_equations!(eqs, domain::WellGroup, model::SimulationModel)
    eqs[:control_equation] = ControlEquationWell()
end

function setup_forces(model::SimulationModel{D}; control = nothing, limits = nothing, set_default_limits = true) where {D <: WellGroup}
    T = Dict{Symbol, Any}
    if isnothing(control)
        control = T()
    end
    wells = model.domain.well_symbols
    for w in wells
        if !haskey(control, w)
            control[w] = DisabledControl()
        end
    end
    # Initialize limits
    if isnothing(limits)
        limits = T()
    end
    for w in wells
        if set_default_limits
            # Set default limits with reasonable values (e.g. minimum rate limit on bhp producers)
            defaults = default_limits(control[w])
            if haskey(limits, w)
                if !isnothing(defaults)
                    if isnothing(limits[w])
                        limits[w] = defaults
                    else
                        limits[w] = merge(defaults, limits[w])
                    end
                end
            else
                limits[w] = defaults
            end
        else
            # Ensure that all limits exist, but set to nothing if not already present
            if !haskey(limits, w)
                limits[w] = nothing
            end
        end
    end
    return (control = control::AbstractDict, limits = limits::AbstractDict,)
end

function convergence_criterion(model, storage, eq::ControlEquationWell, eq_s, r; dt = 1.0, update_report = missing)
    wells = model.domain.well_symbols
    cfg = storage.state.WellGroupConfiguration
    e = abs.(vec(r))
    names = map(w -> name_equation(w, cfg), wells)
    R = (Abs = (errors = e, names = names), )
    return R
end

function name_equation(name, cfg::WellGroupConfiguration)
    ctrl = cfg.operating_controls[name]
    if ctrl isa InjectorControl
        cs = "I"
    elseif ctrl isa ProducerControl
        cs = "P"
    else
        cs = "X"
    end
    t = translate_target_to_symbol(ctrl.target)
    return "$name ($cs) $t"
end
