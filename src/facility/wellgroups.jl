function Jutul.active_entities(wg::WellGroup, ::Jutul.TrivialGlobalMap, e; for_variables = false)
    return 1:count_entities(wg, e)
end

function Jutul.count_entities(wg::WellGroup, ::Wells)
    return length(wg.well_symbols)
end

function get_well_position(d, symbol)
    return findfirst(isequal(symbol), d.well_symbols)
end

function Jutul.associated_entity(::SurfaceTemperature)
    return Wells()
end

function Jutul.associated_entity(::TotalSurfaceMassRate)
    return Wells()
end

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

function Jutul.update_primary_variable!(state, var::SurfacePhaseRates, state_symbol, model, dx, w)
    v = state[state_symbol]
    symbols = model.domain.well_symbols
    cfg = state.WellGroupConfiguration
    # Injectors can only have strictly positive injection rates,
    # producers can only have strictly negative and disabled controls give zero rate.
    abs_max = Jutul.absolute_increment_limit(var)
    rel_max = Jutul.relative_increment_limit(var)
    function do_update!(wcfg, s, v, dx, ctrl)
        return Jutul.update_value(v, dx)
    end
    function do_update!(wcfg, s, v, dx, ctrl::InjectorControl)
        limit_rate = 0.0
        next = Jutul.update_value(v, dx, abs_max, rel_max, limit_rate, nothing)
        return next
    end
    function do_update!(wcfg, s, v, dx, ctrl::ProducerControl)
        # A significant negative rate is the valid producer control
        limit_rate = 0.0
        next = Jutul.update_value(v, dx, abs_max, rel_max, nothing, limit_rate)
        return next
    end
    function do_update!(wcfg, s, v, dx, ctrl::DisabledControl)
        # Set value to zero since we know it is correct.
        return Jutul.update_value(v, -value(v))
    end
    for i in eachindex(symbols)
        s = symbols[i]
        for ph in axes(v, 1)
            v[ph, i] = do_update!(cfg, s, v[ph, i], w*dx[ph, i], operating_control(cfg, s))
        end
    end
    return state
end


rate_weighted(t) = true
rate_weighted(::BottomHolePressureTarget) = false
rate_weighted(::DisabledTarget) = false


target_scaling(::Any) = 1.0
target_scaling(::BottomHolePressureTarget) = 1e5
target_scaling(::SurfaceVolumeTarget) = 1.0/(3600*24.0) # weight by day
target_scaling(::TotalMassRateTarget) = 1.0/(3600*24.0) # weight by day

Jutul.associated_entity(::ControlEquationWell) = Wells()
Jutul.local_discretization(::ControlEquationWell, i) = nothing

function Jutul.prepare_equation_in_entity!(i, eq::ControlEquationWell, eq_s, state, state0, model, dt)
    well = model.domain.well_symbols[i]
    cond = FacilityVariablesForWell(model, state, well, drop_ad = true)
    cfg = state.WellGroupConfiguration
    ctrl = operating_control(cfg, well)
    limits = current_limits(cfg, well)
    apply_well_limits!(cfg, model, state, limits, ctrl, well, cond)
end

function Jutul.update_equation_in_entity!(v, i, state, state0, eq::ControlEquationWell, model, dt, ldisc = local_discretization(eq, i))
    well = model.domain.well_symbols[i]
    cond = FacilityVariablesForWell(model, state, well)
    cfg = state.WellGroupConfiguration
    ctrl = operating_control(cfg, well)
    eq_val = well_control_equation(ctrl, cond, well, model, state)
    v[1] = eq_val + 0*state.BottomHolePressure[i] + 0*state.SurfacePhaseRates[1, i]
end

Jutul.associated_entity(::SurfaceTemperatureEquation) = Wells()
Jutul.local_discretization(::SurfaceTemperatureEquation, i) = nothing
function Jutul.update_equation_in_entity!(v, i, state, state0, eq::SurfaceTemperatureEquation, model, dt, ldisc = local_discretization(eq, i))
    # Set equal to surface temperature. corresponding well temperatures will be
    # subtracted using corss terms
    v[1] = state.SurfaceTemperature[i]
end

Jutul.associated_entity(::BottomHolePressureEquation) = Wells()
Jutul.local_discretization(::BottomHolePressureEquation, i) = nothing
function Jutul.update_equation_in_entity!(v, i, state, state0, eq::BottomHolePressureEquation, model, dt, ldisc = local_discretization(eq, i))
    # Set equal to bhp. corresponding well top cell pressures will be
    # subtracted using corss terms
    v[1] = state.BottomHolePressure[i]
end

Jutul.associated_entity(::SurfacePhaseRatesEquation) = Wells()
Jutul.local_discretization(::SurfacePhaseRatesEquation, i) = nothing
Jutul.number_of_equations_per_entity(fmodel::SimulationModel, eq::SurfacePhaseRatesEquation) = length(get_phases(fmodel.system))

function Jutul.update_equation_in_entity!(v, i, state, state0, eq::SurfacePhaseRatesEquation, model, dt, ldisc = local_discretization(eq, i))
    # Set equal to bhp. corresponding well top cell pressures will be
    # subtracted using corss terms
    for ph in eachindex(v)
        v[ph] = state.SurfacePhaseRates[ph, i]
    end
    return v
end

# Selection of primary variables
function select_primary_variables!(S, system::FacilitySystem, model::FacilityModel)
    ph = get_phases(system.multiphase)
    S[:BottomHolePressure] = BottomHolePressure()
    S[:SurfacePhaseRates] = SurfacePhaseRates(ph)
    S[:TotalSurfaceMassRate] = TotalSurfaceMassRate()
end

function select_equations!(eqs, system::FacilitySystem, model::FacilityModel)
    eqs[:bottom_hole_pressure_equation] = BottomHolePressureEquation()
    eqs[:surface_phase_rates_equation] = SurfacePhaseRatesEquation()
    eqs[:control_equation] = ControlEquationWell()
end

function setup_forces(model::SimulationModel{D}; control = nothing, limits = nothing, set_default_limits = true, check = false) where {D <: WellGroup}
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
        c = control[w]
        is_inj = c isa InjectorControl
        is_prod = c isa ProducerControl
        if is_inj || is_prod
            if is_inj
                sgn = 1.0
            else
                sgn = -1.0
            end
            if check && !isnothing(limits[w])
                for l in [:rate, :orat, :wrat, :grat, :lrat, :resv]
                    if haskey(limits[w], l)
                        val = limits[w][l]
                        if sign(val) != sgn
                            @warn "Well '$w' has $(l) limit with wrong sign ($(val)) for its control type ($(typeof(c)))."
                        end
                    end
                end
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
