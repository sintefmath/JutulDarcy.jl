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
        if !changed && !well_was_disabled && !(new_is_disabled || old_is_disabled)
            # Handle the case where controls may be identical, but a higher type
            # might be incoming (e.g. use of AD)
            changed = changed || value_has_promoted_type(newctrl.target.value, oldctrl.target.value)
        end
        if changed
            # We have a new control. Any previous control change is invalid.
            # Set both operating and requested control to the new one.
            @debug "Well $key switching from $oldctrl to $newctrl"
            # idx = get_well_position(model.domain, key)
            req_ctrls[key] = newctrl
            op_ctrls[key] = newctrl
            cond = FacilityVariablesForWell(model, storage.state, key, drop_ad = true)
            set_facility_values_for_control!(storage.state, model, newctrl, cfg.limits[key], cond)
        end
        pos = get_well_position(model.domain, key)
        # if q_t isa Vector
        #     q_t[pos] = valid_surface_rate_for_control(q_t[pos], newctrl)
        # end
        if changed && !new_is_disabled
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

function valid_surface_rate_for_control(q_t, ::InjectorControl)
    if q_t < MIN_INITIAL_WELL_RATE
        q_t = Jutul.replace_value(q_t, MIN_INITIAL_WELL_RATE)
    end
    return q_t
end

function valid_surface_rate_for_control(q_t, ::ProducerControl)
    if q_t > -MIN_INITIAL_WELL_RATE
        q_t = Jutul.replace_value(q_t, -MIN_INITIAL_WELL_RATE)
    end
    return q_t
end

function valid_surface_rate_for_control(q_t, ::DisabledControl)
    return Jutul.replace_value(q_t, 0.0)
end

function apply_well_limits!(cfg::WellGroupConfiguration, model, state, limits, control, well::Symbol, cond::FacilityVariablesForWell)
    if control isa DisabledControl
        # Disabled wells cannot change active constraint
        println("Well $well is disabled; skipping limit application.")
        return cfg
    end

    if isnothing(limits)
        # No limits to apply
        return cfg
    end
    println("Checking: $well with limits $limits with $cond")
    old_control = control
    control, changed = check_well_limits(limits, cond, control)
    if changed
        @info "Well $well switching control from $(old_control.target) to $(control.target) due to active limit." limits
        cfg.operating_controls[well] = control
        set_facility_values_for_control!(state, model, control, limits, cond)
        # error()
    end
    println("$well operating $(translate_target_to_symbol(control.target)) with value $(control.target.value)")
    return cfg
end

function set_facility_values_for_control!(state, model::FacilityModel, control, limits, cond::FacilityVariablesForWell)
    idx = cond.idx
    @warn "$(cond.name) Setting facility values for control $control at $idx"
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
        gval = wval = 0.0
        oval = w*abs(limits.orat)
    elseif new_target_symbol == :grat
        oval = wval = 0.0
        gval = w*abs(limits.grat)
    elseif new_target_symbol == :wrat
        oval = gval = 0.0
        wval = w*abs(limits.wrat)
    end
    bhps = state.BottomHolePressure
    bhp_limit = get(limits, :bhp, nothing)
    if new_target_symbol == :bhp
        bhps[idx] = replace_value(bhps[idx], bhp_limit)
    elseif !isnothing(bhp_limit)
        bhp = value(bhps[idx])
        is_above = bhp > bhp_limit
        is_below = bhp < bhp_limit
        if (is_injector && is_above) || (!is_injector && is_below)
            @error "Setting bhp" limits.bhp*(1.0 - 0.01*sgn) is_injector bhp bhp_limit
            bhps[idx] = replace_value(bhps[idx], bhp_limit*(1.0 - 0.01*sgn))
        end
    end
    ev = 1e-8
    phase_rates = state.SurfacePhaseRates
    function set_phase_rate!(ph, val)
        phase_rates[ph, idx] = replace_value(phase_rates[ph, idx], sgn*val)
    end

    q_t = state.TotalSurfaceMassRate
    q_t[idx] = valid_surface_rate_for_control(q_t[idx], control)

    do_print = cond.name == :PRODU21
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
    if do_print
        @info "??" value(phase_rates[:, idx]) value(bhps[idx]) cond
    end
    return state
end

function check_well_limits(limits, cond, control)
    current_target = control.target
    next_target = current_target
    current_name = translate_target_to_symbol(current_target)
    if control isa InjectorControl
        for (name, limit_value) in pairs(limits)
            if name == current_name
                continue
            end

            if name == :bhp
                # Injector BHP limit is an upper limit, and pressure is positive
                if cond.bottom_hole_pressure > limit_value
                    @error "INJECTOR BHP TOO HIGH" cond.bottom_hole_pressure limit_value
                    next_target = BottomHolePressureTarget(limit_value)
                    break
                end
            else
                # Injector limits are upper limits on rates
                if name == :rate
                    rate = cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate
                    if rate > limit_value
                        next_target = TotalRateTarget(limit_value)
                        break
                    end
                elseif name == :wrat
                    wrat = cond.surface_aqueous_rate
                    if wrat > limit_value
                        next_target = SurfaceWaterRateTarget(limit_value)
                        break
                    end
                elseif name == :orat
                    orat = cond.surface_liquid_rate
                    if orat > limit_value
                        next_target = SurfaceOilRateTarget(limit_value)
                        break
                    end
                elseif name == :lrat
                    lrat = cond.surface_aqueous_rate + cond.surface_liquid_rate
                    if lrat > limit_value
                        next_target = SurfaceLiquidRateTarget(limit_value)
                        break
                    end
                elseif name == :grat
                    grat = cond.surface_vapor_rate
                    if grat > limit_value
                        next_target = SurfaceGasRateTarget(limit_value)
                        break
                    end
                elseif name == :rate_lower
                    rate = cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate
                    if rate < limit_value
                        next_target = TotalRateTarget(limit_value)
                        break
                    end
                else
                    error("Unsupported well producer constraint/limit $k")
                end
            end
        end
    else
        control::ProducerControl
        for (name, limit_value) in pairs(limits)
            if name == current_name
                continue
            end
            if name == :bhp
                # Producer BHP limit is a lower limit, and pressure is positive
                if cond.bottom_hole_pressure < limit_value
                    @error "PRODUCER BHP TOO LOW" cond.bottom_hole_pressure limit_value
                    next_target = BottomHolePressureTarget(limit_value)
                    break
                end
            else
                # Producer rates are negative by convention. Producer rate
                # limits are upper limits.
                limit_value = abs(limit_value)
                if name == :lrat
                    lrat = -(cond.surface_aqueous_rate + cond.surface_liquid_rate)
                    if lrat > limit_value
                        next_target = SurfaceLiquidRateTarget(-limit_value)
                        break
                    end
                elseif name == :orat
                    orat = -cond.surface_liquid_rate
                    if orat > limit_value
                        next_target = SurfaceOilRateTarget(-limit_value)
                        break
                    end
                elseif name == :wrat
                    wrat = -cond.surface_aqueous_rate
                    if wrat > limit_value
                        next_target = SurfaceWaterRateTarget(-limit_value)
                        break
                    end
                elseif name == :grat
                    grat = -cond.surface_vapor_rate
                    if grat > limit_value
                        next_target = SurfaceGasRateTarget(-limit_value)
                        break
                    end
                elseif name == :rate
                    rate = -(cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate)
                    if rate > limit_value
                        next_target = TotalRateTarget(-limit_value)
                        break
                    end
                elseif name == :rate_lower
                    rate = -(cond.surface_aqueous_rate + cond.surface_liquid_rate + cond.surface_vapor_rate)
                    if rate < limit_value
                        next_target = TotalRateTarget(-limit_value)
                        break
                    end
                else
                    error("Unsupported well injector constraint/limit $k")
                end
            end
        end
    end
    changed = next_target != current_target
    if changed
        @error "Well changed" next_target current_target cond
        control = replace_target(control, next_target)
    end
    return (control, next_target != current_target)
end


# function apply_well_limit!(cfg::WellGroupConfiguration, target, wmodel, wstate, well::Symbol, density_s, volume_fraction_s, total_mass_rate, current_lims = current_limits(cfg, well))
#     if !isnothing(current_lims)
#         ctrl = operating_control(cfg, well)
#         @tic "limits" target, changed, current_val, limit_val, lim_type = check_active_limits(ctrl, target, current_lims, wmodel, wstate, well, density_s, volume_fraction_s, total_mass_rate)
#         if changed
#             old = cfg.operating_controls[well].target
#             next_control = replace_target(ctrl, target)
#             cfg.operating_controls[well] = next_control
#             if lim_type == :lower
#                 lb = "is below"
#             else
#                 lb = "is above"
#             end
#             limstr = "Current value $(current_val) $lb $lim_type limit $(limit_val)."
#             header = "$well is switching control due to $lim_type limit for $(typeof(target)) of $limit_val"
#             @debug "$(header)\n$(limstr)\nOld value: $old\nNew value: $target." next_control
#         end
#     end
#     return target
# end

# function check_active_limits(control, target, limits, wmodel, wstate, well::Symbol, density_s, volume_fraction_s, total_mass_rate)
#     changed = false
#     cval = tval = NaN
#     is_lower = false
#     for (name, val) in pairs(limits)
#         if isfinite(first(val))
#             (target_limit, is_lower) = translate_limit(control, name, val)
#             ok, cval, tval = check_limit(control, target_limit, target, is_lower, total_mass_rate, wmodel, wstate, density_s, volume_fraction_s)
#             if !ok
#                 changed = true
#                 target = target_limit
#                 break
#             end
#         end
#     end
#     if is_lower
#         lim_type = :lower
#     else
#         lim_type = :upper
#     end
#     return (target, changed, cval, tval, lim_type)
# end

function well_control_equation(ctrl::DisabledControl, cond, well, model, state)
    target = ctrl.target
    target::DisabledTarget
    # Equation is just the sum of absolute rates to force them all to zero.
    return well_target_value(ctrl, target, cond, well, model, state)
end

function well_control_equation(ctrl, cond, well, model, state)
    target = ctrl.target
    val_t = target.value
    val = well_target_value(ctrl, target, cond, well, model, state)
    scale = target_scaling(target)
    # @info "Target: $(val_t) - $(value(val)) for $cond and $(ctrl.target)"
    return (val - val_t)/scale
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

# """
#     translate_limit(control::ProducerControl, name, val)

# Translates the limit for a given control parameter in a `ProducerControl`.

# # Arguments
# - `control::ProducerControl`: The control object containing the parameters to be translated into limit.
# - `name`: The name of the parameter whose limit is to be translated.
#     - `:bhp`: Bottom hole pressure.
#     - `:orat`: Surface oil rate.
#     - `:lrat`: Surface liquid (water + oil) rate.
#     - `:grat`: Surface gas rate.
#     - `:wrat`: Surface water rate.
#     - `:rate`: Total volumetric surface rate (upper limit).
#     - `:rate_upper`: Total volumetric surface rate (upper limit).
#     - `:rate_lower`: Total volumetric surface rate (lower limit).
#     - `:resv`: Reservoir voidage.
#     - `:mrat`: Total mass rate (upper limit).
# - `val`: The value to which the limit is to be translated.

# # Returns
# - The translated limit value for the specified control parameter.
# """
# function translate_limit(control::ProducerControl, name, val)
#     # Note: Negative sign convention for production.
#     # A lower absolute bound on a rate
#     # |q| > |lim| -> q < lim if both sides are negative
#     # means that we specify is_lower for upper limits and the other
#     # way around for lower limits, when dealing with rates.
#     is_lower = true
#     if name == :bhp
#         # Upper limit, pressure
#         target_limit = BottomHolePressureTarget(val)
#         # Pressures are positive, this is a true lower bound
#         is_lower = true
#     elseif name == :orat
#         # Upper limit, surface oil rate
#         target_limit = SurfaceOilRateTarget(val)
#     elseif name == :lrat
#         # Upper limit, surface liquid (water + oil) rate
#         target_limit = SurfaceLiquidRateTarget(val)
#     elseif name == :grat
#         # Upper limit, surface gas rate
#         target_limit = SurfaceGasRateTarget(val)
#     elseif name == :wrat
#         # Upper limit, surface water rate
#         target_limit = SurfaceWaterRateTarget(val)
#     elseif name == :rate || name == :rate_upper
#         # Upper limit, total volumetric surface rate
#         target_limit = TotalRateTarget(val)
#     elseif name == :rvolrat
#         target_limit = ReservoirVolumeRateTarget(val)
#     elseif name == :rate_lower
#         # Lower limit, total volumetric surface rate. This is useful
#         # disabling producers if they would otherwise start to inject.
#         target_limit = TotalRateTarget(val)
#         is_lower = false
#     elseif name == :resv
#         v, w = val
#         target_limit = ReservoirVoidageTarget(v, w)
#     elseif name == :mrat
#         # Upper limit, total mass rate
#         target_limit = TotalMassRateTarget(val)
#     else
#         error("$name limit not supported for well acting as producer.")
#     end
#     return (target_limit, is_lower)
# end

# """
#     translate_limit(control::InjectorControl, name, val)

# Translate the limit for a given `InjectorControl` object.

# # Arguments
# - `control::InjectorControl`: The control object for which the limit is being translated.
# - `name`: The name of the limit to be translated.
#     - `:bhp`: Bottom hole pressure.
#     - `:rate`: Total volumetric surface rate (upper limit).
#     - `:rate_upper`: Total volumetric surface rate (upper limit).
#     - `:rate_lower`: Total volumetric surface rate (lower limit).
#     - `:resv_rate`: Total volumetric reservoir rate.
#     - `:mrat: Total mass rate(upper limit)
# - `val`: The value associated with the limit.

# # Returns
# - The translated limit value.
# """
# function translate_limit(control::InjectorControl, name, val)
#     is_lower = false
#     if name == :bhp
#         # Upper limit, pressure
#         target_limit = BottomHolePressureTarget(val)
#     elseif name == :rate || name == :rate_upper || name == :wrat || name == :orat || name == :lrat || name == :grat
#         # Upper limit, total volumetric surface rate
#         target_limit = TotalRateTarget(val)
#     elseif name == :rate_lower
#         # Lower limit, total volumetric surface rate
#         target_limit = TotalRateTarget(val)
#         is_lower = true
#     elseif name == :resv_rate
#         # Upper limit, total volumetric reservoir rate
#         target_limit = TotalReservoirRateTarget(val)
#     elseif name == :mrat
#         # Upper limit, total mass rate
#         target_limit = TotalMassRateTarget(val)
#     else
#         error("$name limit not supported for well acting as injector.")
#     end
#     return (target_limit, is_lower)
# end

# function check_limit(current_control, target_limit, target, is_lower::Bool, q_t, source_model, well_state, rhoS, S)
#     if typeof(target_limit) == typeof(target)
#         # We are already operating at this target and there is no need to check.
#         ok = true
#         current_val = limit_val = NaN
#     else
#         current_val = value(well_target_value(q_t, current_control, target_limit, source_model, well_state, rhoS, S))
#         limit_val = target_limit.value
#         ϵ = 1e-6
#         if is_lower
#             # Limit is lower bound, check that we are above...
#             ok = current_val >= (1 + ϵ)*limit_val
#         else
#             ok = current_val <= (1 - ϵ)*limit_val
#         end
#     end
#     return (ok, current_val, limit_val)
# end


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

# """
# Well target contribution from well itself (disabled, zero value)
# """
# function well_target(control, target::DisabledTarget, well_model, well_state, rhoS, S)
#     return 0.0
# end

# """
# Well target contribution from well itself (bhp)
# """
# function well_target(control, target::BottomHolePressureTarget, well_model, well_state, rhoS, S)
#     return bottom_hole_pressure(well_state)
# end

# """
# Well target contribution from well itself (surface volume, injector)
# """
# function well_target(control::InjectorControl, target::SurfaceVolumeTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     t_phases = lumped_phases(target)
#     w_phases = get_phases(well_model.system)
#     t = 0.0
#     for (ix, mix) in control.phases
#         if w_phases[ix] in t_phases
#             t += mix
#         end
#     end
#     return t/control.mixture_density
# end


# """
# Well target contribution from well itself (reservoir volume, injector)
# """
# function well_target(control::InjectorControl, target::TotalReservoirRateTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     w_phases = get_phases(well_model.system)
#     t = 0.0
#     rho = well_state.PhaseMassDensities
#     for (ix, mix) in control.phases
#         t += mix/rho[ix, 1]
#     end
#     return t
# end

# """
# Well target contribution from well itself (surface volume, injector)
# """
# function well_target(control::InjectorControl, target::TotalRateTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     if abs(target.value) == MIN_ACTIVE_WELL_RATE
#         # Special meaning: Limits set at the CONST minimum value are really absolute mass limits.
#         w = 1.0
#     else
#         w = 1.0/control.mixture_density
#     end
#     return w
# end

# """
# Well target contribution from well itself (surface volume, producer)
# """
# function well_target(control::ProducerControl, target::SurfaceVolumeTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     phases = get_phases(well_model.system)
#     Tw = eltype(surface_volume_fractions)
#     if abs(target.value) == MIN_ACTIVE_WELL_RATE
#         # Special meaning: Limits set at the CONST minimum value are really absolute mass limits.
#         w = 1.0
#     else
#         # Compute total density at surface conditions by weighting phase volumes at surf
#         ρ_tot = zero(Tw)
#         for (ρ, V) in zip(surface_densities, surface_volume_fractions)
#             ρ_tot += ρ*V
#         end
#         # Divide by total density to get total volume at surface, then multiply that by surface volume fraction
#         w = zero(Tw)
#         if isa(target, TotalRateTarget)
#             for i in eachindex(phases)
#                 @inbounds V = surface_volume_fractions[i]
#                 w += V
#             end
#         else
#             lp = lumped_phases(target)
#             for (i, ph) in enumerate(phases)
#                 if ph in lp
#                     @inbounds V = surface_volume_fractions[i]
#                     w += V
#                 end
#             end
#         end
#         w = w/ρ_tot
#     end
#     return w
# end


# """
# Well target contribution for a producer with TotalMassRateTarget.
# Always returns 1.0 so that `t = q_t` in `target_actual_pair`, i.e., the residual
# directly compares the actual total mass flow (kg/s) with the target value.
# """

# function well_target(control::ProducerControl, target::TotalMassRateTarget,
#                      well_model, well_state, surface_densities, surface_volume_fractions)
#     return 1.0
# end


# """
# Well target contribution from injector with TotalMassRateTarget.

# This target enforces a specified total mass injection rate (kg/s) for the well.

# # Example usage:
# rate_target = TotalMassRateTarget(1.0)

# # Create injector control:
# I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = 630.0)
# # or, mass rate only (density defaults to 1.0 kg/m^3;
# I_ctrl = InjectorControl(rate_target, [0.0, 1.0])
# """

# function well_target(control::InjectorControl, target::TotalMassRateTarget,
#                      well_model, well_state, surface_densities, surface_volume_fractions)                   
#     return 1.0
# end


# """
# Well target contribution from well itself (RESV, producer)
# """
# function well_target(control::ProducerControl, target::ReservoirVoidageTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     Tw = eltype(surface_volume_fractions)
#     ρ_tot = zero(Tw)
#     for (ρ, V) in zip(surface_densities, surface_volume_fractions)
#         ρ_tot += ρ*V
#     end
#     w = zero(Tw)
#     for (i, S) in enumerate(surface_volume_fractions)
#         w += S*target.weights[i]
#     end
#     return w/ρ_tot
# end

# function well_target(control::ProducerControl, target::ReservoirVolumeRateTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     rho = well_state.PhaseMassDensities
#     s = well_state.Saturations
#     total_density = zero(eltype(rho))
#     wc = JutulDarcy.well_top_node()
#     for ph in axes(rho, 1)
#         total_density += rho[ph, wc]*s[ph, wc]
#     end
#     if haskey(well_state, :TotalSaturation)
#         # sT = well_state.TotalSaturation[wc]
#         # total_density *= sT
#     end
#     return 1.0/total_density
# end

# function well_target(control::InjectorControl, target::ReinjectionTarget, well_model, well_state, surface_densities, surface_volume_fractions)
#     w = 1.0/control.mixture_density
#     return w
# end

# function well_target_value(q_t, control, target, source_model, well_state, rhoS, S)
#     v = well_target(control, target, source_model, well_state, rhoS, S)
#     if rate_weighted(target)
#         v *= value(q_t)
#     end
#     return v
# end

