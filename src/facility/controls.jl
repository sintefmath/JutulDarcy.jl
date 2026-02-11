include("special_controls/special_controls.jl")

function Jutul.initialize_extra_state_fields!(state, domain::WellGroup, model; T = Float64)
    # Insert structure that holds well control (limits etc) that is then updated before each step
    state[:WellGroupConfiguration] = WellGroupConfiguration(domain.well_symbols)
end

function Jutul.update_before_step_multimodel!(storage_g, model_g::MultiModel, model::WellGroupModel, dt, forces_g, key;
        time = NaN,
        recorder = ProgressRecorder(),
        update_explicit = true
    )
    function value_has_promoted_type(newval::T1, oldval::T2) where {T1, T2}
        return T1 != T2 && promote_type(T1, T2) == T1
    end
    forces = forces_g[key]
    # Set control to whatever is on the forces
    storage = storage_g[key]
    cfg = storage.state.WellGroupConfiguration
    q_t = storage.state.TotalSurfaceMassRate
    op_ctrls = cfg.operating_controls
    req_ctrls = cfg.requested_controls
    # Set limits
    for key in keys(forces.limits)
        cfg.limits[key] = forces.limits[key]
    end
    # Set group controls if present in forces
    if haskey(forces, :group_controls)
        for (gname, gctrl) in forces.group_controls
            cfg.group_controls[gname] = gctrl
        end
    end
    current_step = recorder.recorder.step
    # Set operational controls
    for key in keys(forces.control)
        # If the requested control in forces differ from the one we are presently using, we need to switch.
        # Otherwise, stay the course.
        rmodel = model_g[:Reservoir]
        rstate = storage_g.Reservoir.state
        newctrl, changed = realize_control_for_reservoir(rstate, forces.control[key], rmodel, dt)
        oldctrl = req_ctrls[key]
        is_new_step = cfg.step_index != current_step
        disabled = DisabledControl()
        new_is_disabled = newctrl isa DisabledControl
        old_is_disabled = oldctrl isa DisabledControl
        well_was_disabled = op_ctrls[key] == disabled && !new_is_disabled
        changed = (is_new_step && newctrl != oldctrl) || well_was_disabled
        has_group_target = newctrl.target isa GroupTarget
        any_disabled = new_is_disabled || old_is_disabled
        if !changed && !well_was_disabled && !any_disabled && !has_group_target
            # Handle the case where controls may be identical, but a higher type
            # might be incoming (e.g. use of AD)
            changed = changed || value_has_promoted_type(newctrl.target.value, oldctrl.target.value)
        end
        if changed
            # We have a new control. Any previous control change is invalid.
            # Set both operating and requested control to the new one.
            @debug "Well $key switching from $oldctrl to $newctrl"
            req_ctrls[key] = newctrl
            op_ctrls[key] = newctrl
        end
        pos = get_well_position(model.domain, key)
        if q_t isa Vector
            q_t[pos] = valid_surface_rate_for_control(q_t[pos], newctrl)
        end
        if changed && !new_is_disabled && !has_group_target
            if isnothing(cfg.limits[key])
                cfg.limits[key] = as_limit(newctrl.target)
            else
                cfg.limits[key] = merge(cfg.limits[key], as_limit(newctrl.target))
            end
        end
    end
    cfg.step_index = current_step
    for wname in model.domain.well_symbols
        wmodel = model_g[wname]
        wstate = storage_g[wname].state
        forces_w = forces_g[wname]
        if isnothing(forces_w) || !haskey(forces_w, :mask)
            mask = nothing
        else
            mask = forces_w.mask
        end
        rmodel = model_g[:Reservoir]
        rstate = storage_g.Reservoir.state
        update_before_step_well!(wstate, wmodel, rstate, rmodel, op_ctrls[wname], mask, update_explicit = update_explicit)
    end
end


function valid_surface_rate_for_control(q_t, ::InjectorControl, val = MIN_INITIAL_WELL_RATE)
    if q_t < val
        q_t = Jutul.replace_value(q_t, val)
    end
    return q_t
end

function valid_surface_rate_for_control(q_t, ::ProducerControl, val = MIN_INITIAL_WELL_RATE)
    if q_t > -val
        q_t = Jutul.replace_value(q_t, -val)
    end
    return q_t
end

function valid_surface_rate_for_control(q_t, ::DisabledControl)
    return Jutul.replace_value(q_t, 0.0)
end

function apply_well_limits!(cfg::WellGroupConfiguration, model, state, limits, control, well::Symbol, cond::FacilityVariablesForWell)
    if control isa DisabledControl
        # Disabled wells cannot change active constraint
        return cfg
    end

    if isnothing(limits)
        # No limits to apply
        return cfg
    end
    old_control = control
    control, changed = check_well_limits(limits, cond, control)
    if changed
        @debug "Well $well switching control from $(old_control.target) to $(control.target) due to active limit." limits cond
        cfg.operating_controls[well] = control
    end
    return cfg
end

function set_facility_values_for_control!(state, model::FacilityModel, control, limits, cond::FacilityVariablesForWell)
    idx = cond.idx
    q_t = state.TotalSurfaceMassRate
    q_t[idx] = valid_surface_rate_for_control(q_t[idx], control)
    return
    # @warn "$(cond.name) Setting facility values for control $control at $idx"
    new_target_symbol = translate_target_to_symbol(control.target)
    is_injector = control isa InjectorControl
    if is_injector
        sgn = 1.0
    else
        sgn = -1.0
    end
    sys = model.system.multiphase
    nph = number_of_phases(sys)
    rhos = reference_densities(sys)

    w_idx = phase_index(sys, AqueousPhase())
    o_idx = phase_index(sys, LiquidPhase())
    g_idx = phase_index(sys, VaporPhase())

    has_water = !isnothing(w_idx)
    has_oil = !isnothing(o_idx)
    has_gas = !isnothing(g_idx)

    wval = abs(cond.surface_aqueous_rate)
    oval = abs(cond.surface_liquid_rate)
    gval = abs(cond.surface_vapor_rate)

    # Upper limits
    w = 0.99
    w = 1.0
    if new_target_symbol == :lrat
        lrat_limit = w*abs(limits.lrat)
        if has_water && has_oil
            wval = lrat_limit/2.0
            oval = lrat_limit/2.0
        elseif has_water
            wval = lrat_limit
        elseif has_oil
            oval = lrat_limit
        else
            error("No water or oil phase present to apply liquid rate limit.")
        end
        gval = 0.0
    elseif new_target_symbol == :rate
        # Injector should be different!
        rate = w*abs(limits.rate)
        if is_injector
            wval = oval = gval = 0.0
            for (ph, mix) in control.phases
                if ph == w_idx
                    wval = mix*rate
                elseif ph == o_idx
                    oval = mix*rate
                elseif ph == g_idx
                    gval = mix*rate
                end
            end
        else
            if has_water
                wval = rhos[w_idx]*rate/nph
            else
                wval = 0.0
            end
            if has_oil
                oval = rhos[o_idx]*rate/nph
            else
                oval = 0.0
            end
            if has_gas
                gval = rhos[g_idx]*rate/nph
            end
            total = abs(wval + oval + gval)
            wval = wval/total
            oval = oval/total
            gval = gval/total
        end
    elseif new_target_symbol == :orat
        # gval = wval = 0.0
        oval = w*abs(limits.orat)
    elseif new_target_symbol == :grat
        # oval = wval = 0.0
        gval = w*abs(limits.grat)
    elseif new_target_symbol == :wrat
        # oval = gval = 0.0
        wval = w*abs(limits.wrat)
    end
    bhps = state.BottomHolePressure
    bhp_limit = get(limits, :bhp, nothing)
    if new_target_symbol == :bhp
        # bhps[idx] = replace_value(bhps[idx], bhp_limit)
    elseif !isnothing(bhp_limit)
        bhp = value(bhps[idx])
        is_above = bhp > bhp_limit
        is_below = bhp < bhp_limit
        if (is_injector && is_above) || (!is_injector && is_below)
            # @error "Setting bhp" limits.bhp*(1.0 - 0.01*sgn) is_injector bhp bhp_limit
            # bhps[idx] = replace_value(bhps[idx], bhp_limit*(1.0 - 0.01*sgn))
        end
    end
    ev = 1e-8
    phase_rates = state.SurfacePhaseRates
    function set_phase_rate!(ph, val)
        # phase_rates[ph, idx] = replace_value(phase_rates[ph, idx], sgn*val)
    end


    do_print = cond.name == :PRODU21 && false
    if do_print
        @info "Setting values for $new_target_symbol $(limits[new_target_symbol])" wval oval gval value(bhps[idx]) limits
    end
    if has_water
        set_phase_rate!(w_idx, max(wval - ev, 0))
    end
    if has_oil
        set_phase_rate!(o_idx, max(oval - ev, 0))
    end
    if has_gas
        set_phase_rate!(g_idx, max(gval - ev, 0))
    end
    new_mass_rate = 0.0
    rhos = reference_densities(sys)
    for ph in 1:number_of_phases(sys)
        new_mass_rate += rhos[ph]*value(phase_rates[ph, idx])
    end
    qt0 = value(q_t[idx])
    # q_t[idx] = replace_value(q_t[idx], new_mass_rate)
    if do_print
        @info "???" new_mass_rate qt0
        @info "??" value(phase_rates[:, idx]) value(bhps[idx]) cond
    end
    return state
end

function check_well_limits(limits, cond, control)
    # current_target = control.target
    # next_target = current_target
    current_name = translate_target_to_symbol(control.target)
    next_target = missing
    changed = false
    if haskey(limits, current_name)
        control_check = check_well_limit(current_name, limits[current_name], cond, control)
        control_ok = ismissing(control_check)
    else
        control_ok = true
    end
    if control_ok
        for (name, limit_value) in pairs(limits)
            if name == current_name
                continue
            end
            next_target = check_well_limit(name, limit_value, cond, control)
            if !ismissing(next_target)
                changed = true
                break
            end
        end
    end
    if changed
        # @error "Well changed" next_target current_target cond
        control = replace_target(control, next_target)
    end
    return (control, changed)
end

"""
    check_well_limit(name::Symbol, limit_value, cond, control::InjectorControl)

Check constrants/limits for an injector well control. If a limit is violated, return
a new target to switch to. If no limits are violated, return `missing`.
"""
function check_well_limit(name::Symbol, limit_value, cond, control::InjectorControl)
    next_target = missing
    if name == :bhp
        # Injector BHP limit is an upper limit, and pressure is positive
        if cond.bottom_hole_pressure > limit_value
            # @error "INJECTOR BHP TOO HIGH" cond.bottom_hole_pressure limit_value
            next_target = BottomHolePressureTarget(limit_value)
        end
    else
        # Injector limits are upper limits on rates
        if name == :rate
            rate = cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate
            if rate > limit_value
                next_target = TotalRateTarget(limit_value)
            end
        elseif name == :wrat
            wrat = cond.surface_aqueous_rate
            if wrat > limit_value
                next_target = SurfaceWaterRateTarget(limit_value)
            end
        elseif name == :orat
            orat = cond.surface_liquid_rate
            if orat > limit_value
                next_target = SurfaceOilRateTarget(limit_value)
            end
        elseif name == :lrat
            lrat = cond.surface_aqueous_rate + cond.surface_liquid_rate
            if lrat > limit_value
                next_target = SurfaceLiquidRateTarget(limit_value)
            end
        elseif name == :grat
            grat = cond.surface_vapor_rate
            if grat > limit_value
                next_target = SurfaceGasRateTarget(limit_value)
            end
        elseif name == :rate_lower
            rate = cond.total_mass_rate
            if rate < limit_value
                next_target = TotalRateTarget(limit_value)
            end
        elseif name == :reinjection
            rate = cond.total_mass_rate
            if rate > limit_value
                next_target = TotalMassRateTarget(limit_value)
            end
        else
            error("Unsupported well producer constraint/limit $name")
        end
    end
    return next_target
end

"""
    check_well_limit(name::Symbol, limit_value, cond::FacilityVariablesForWell, control::ProducerControl)

Check constrants/limits for a producer well control. If a limit is violated, return
a new target to switch to. If no limits are violated, return `missing`.
"""
function check_well_limit(name::Symbol, limit_value, cond::FacilityVariablesForWell, control::ProducerControl)
    next_target = missing
    if name == :bhp
        # Producer BHP limit is a lower limit, and pressure is positive
        if cond.bottom_hole_pressure < limit_value
            next_target = BottomHolePressureTarget(limit_value)
        end
    else
        # Producer rates are negative by convention. Producer rate
        # limits are upper limits.
        if name == :resv
            state_avg = limit_value[2]
            limit_value = abs(limit_value[1])
            q_w = -cond.surface_aqueous_rate
            q_o = -cond.surface_liquid_rate
            q_g = -cond.surface_vapor_rate
            resv = compute_total_resv_rate(state_avg; qw = q_w, qg = q_g, qo = q_o)
            if resv > limit_value
                next_target = ReservoirVoidageTarget(-limit_value, state_avg)
            end
        else
            limit_value = abs(limit_value)
            if name == :lrat
                lrat = -(cond.surface_aqueous_rate + cond.surface_liquid_rate)
                if lrat > limit_value
                    next_target = SurfaceLiquidRateTarget(-limit_value)
                end
            elseif name == :orat
                orat = -cond.surface_liquid_rate
                if orat > limit_value
                    next_target = SurfaceOilRateTarget(-limit_value)
                end
            elseif name == :wrat
                wrat = -cond.surface_aqueous_rate
                if wrat > limit_value
                    next_target = SurfaceWaterRateTarget(-limit_value)
                end
            elseif name == :grat
                grat = -cond.surface_vapor_rate
                if grat > limit_value
                    next_target = SurfaceGasRateTarget(-limit_value)
                end
            elseif name == :rate
                rate = -(cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate)
                if rate > limit_value
                    next_target = TotalRateTarget(-limit_value)
                end
            elseif name == :rate_lower
                rate = -cond.total_mass_rate
                if rate < limit_value
                    next_target = TotalRateTarget(-limit_value)
                end
            else
                error("Unsupported well injector constraint/limit $name")
            end
        end
    end
    return next_target
end

function well_control_equation(ctrl::DisabledControl, cond::FacilityVariablesForWell, well::Symbol, model, state)
    target = ctrl.target
    target::DisabledTarget
    # Equation is just the sum of absolute rates to force them all to zero.
    return well_target_value(ctrl, target, cond, well, model, state)
end

"""
    well_control_equation(ctrl, cond, well, model, state)

Compute the control equation for the given well control `ctrl`. This function
returns the difference between the current value of the target and
the target value itself, scaled by a target scaling factor.
"""
function well_control_equation(ctrl, cond, well, model, state)
    target = ctrl.target
    if target isa GroupTarget
        return well_control_equation_group(ctrl, target, cond, well, model, state)
    end
    val_t = get_control_target_value(target, model, state)
    val = well_target_value(ctrl, target, cond, well, model, state)
    scale = target_scaling(target)
    return (val - val_t)/scale
end

"""
    well_control_equation_group(ctrl, group_target, cond, well, model, state)

Compute the control equation for a well under group control. The well's actual
target type and value are looked up from the group's control stored in
`WellGroupConfiguration`. For rate-weighted group targets, each well targets
its allocation factor times the group target value. For non-rate targets
(e.g. BHP), the allocation factor is not applied.
"""
function well_control_equation_group(ctrl, group_target::GroupTarget, cond::FacilityVariablesForWell, well::Symbol, model, state)
    group_name = group_target.group
    cfg = state.WellGroupConfiguration
    group_ctrl = cfg.group_controls[group_name]
    actual_target = group_ctrl.target
    val_t = get_control_target_value(actual_target, model, state)
    if rate_weighted(actual_target)
        groups = model.domain.groups
        group_wells = groups[group_name]
        nw = length(group_wells)
        af = isnothing(group_target.allocation_factor) ? 1.0/nw : group_target.allocation_factor
        val_t = val_t * af
    end
    val = well_target_value(ctrl, actual_target, cond, well, model, state)
    scale = target_scaling(actual_target)
    return (val - val_t)/scale
end

function get_control_target_value(target::WellTarget, model, state)
    return target.value
end

function get_control_target_value(target::ReinjectionTarget, model, state)
    qtot = 0.0
    for well in target.wells
        cond_well = FacilityVariablesForWell(model, state, well)
        qtot -= cond_well.total_mass_rate
    end
    return max(qtot, MIN_ACTIVE_WELL_RATE)
end

"""
    well_target_value(ctrl::WellControlForce, target::WellTarget, cond::FacilityVariablesForWell, well::Symbol, model::Facility_model, facility_state)

Compute the value of the well target for the given control and facility
conditions. For example, for the `BottomHolePressureTarget`, this function
returns the bottom hole pressure of the well as given in `cond`. This value is
then used in the control equation to compute the difference between the current
and target values.
"""
function well_target_value

end

function well_target_value(ctrl::DisabledControl, target::DisabledTarget, cond, well, model, state)
    return abs(cond.total_mass_rate) + abs(cond.surface_aqueous_rate) + abs(cond.surface_liquid_rate) + abs(cond.surface_vapor_rate)
end

function well_target_value(ctrl, target::TotalRateTarget, cond, well, model, state)
    rate = cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate
    return rate
end

function well_target_value(ctrl, target::BottomHolePressureTarget, cond, well, model, state)
    return cond.bottom_hole_pressure
end

function well_target_value(ctrl, target::SurfaceOilRateTarget, cond, well, model, state)
    return cond.surface_liquid_rate
end

function well_target_value(ctrl, target::SurfaceWaterRateTarget, cond, well, model, state)
    return cond.surface_aqueous_rate
end

function well_target_value(ctrl, target::SurfaceGasRateTarget, cond, well, model, state)
    return cond.surface_vapor_rate
end

function well_target_value(ctrl, target::SurfaceLiquidRateTarget, cond, well, model, state)
    return cond.surface_liquid_rate + cond.surface_aqueous_rate
end

function well_target_value(ctrl, target::TotalMassRateTarget, cond, well, model, state)
    return cond.total_mass_rate
end

function well_target_value(ctrl, target::ReinjectionTarget, cond, well, model, state)
    return cond.total_mass_rate
end

function facility_surface_mass_rate_for_well(model::SimulationModel, wsym, fstate; effective::Bool = false)
    pos = get_well_position(model.domain, wsym)
    q_t = fstate.TotalSurfaceMassRate[pos]
    if effective
        control = fstate.WellGroupConfiguration.operating_controls[wsym]
        q_t = effective_surface_rate(q_t, control)
    end
    return q_t
end

bottom_hole_pressure(ws) = ws.Pressure[1]

