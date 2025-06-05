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
            req_ctrls[key] = newctrl
            op_ctrls[key] = newctrl
        end
        pos = get_well_position(model.domain, key)
        if q_t isa Vector
            q_t[pos] = valid_surface_rate_for_control(q_t[pos], newctrl)
        end
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

function apply_well_limit!(cfg::WellGroupConfiguration, target, wmodel, wstate, well::Symbol, density_s, volume_fraction_s, total_mass_rate, current_lims = current_limits(cfg, well))
    if !isnothing(current_lims)
        ctrl = operating_control(cfg, well)
        @tic "limits" target, changed, current_val, limit_val, lim_type = check_active_limits(ctrl, target, current_lims, wmodel, wstate, well, density_s, volume_fraction_s, total_mass_rate)
        if changed
            old = cfg.operating_controls[well].target
            next_control = replace_target(ctrl, target)
            cfg.operating_controls[well] = next_control
            if lim_type == :lower
                lb = "is below"
            else
                lb = "is above"
            end
            limstr = "Current value $(current_val) $lb $lim_type limit $(limit_val)."
            header = "$well is switching control due to $lim_type limit for $(typeof(target)) of $limit_val"
            @debug "$(header)\n$(limstr)\nOld value: $old\nNew value: $target." next_control
        end
    end
    return target
end

function check_active_limits(control, target, limits, wmodel, wstate, well::Symbol, density_s, volume_fraction_s, total_mass_rate)
    changed = false
    cval = tval = NaN
    is_lower = false
    for (name, val) in pairs(limits)
        if isfinite(first(val))
            (target_limit, is_lower) = translate_limit(control, name, val)
            ok, cval, tval = check_limit(control, target_limit, target, is_lower, total_mass_rate, wmodel, wstate, density_s, volume_fraction_s)
            if !ok
                changed = true
                target = target_limit
                break
            end
        end
    end
    if is_lower
        lim_type = :lower
    else
        lim_type = :upper
    end
    return (target, changed, cval, tval, lim_type)
end

"""
    translate_limit(control::ProducerControl, name, val)

Translates the limit for a given control parameter in a `ProducerControl`.

# Arguments
- `control::ProducerControl`: The control object containing the parameters to be translated into limit.
- `name`: The name of the parameter whose limit is to be translated.
    - `:bhp`: Bottom hole pressure.
    - `:orat`: Surface oil rate.
    - `:lrat`: Surface liquid (water + oil) rate.
    - `:grat`: Surface gas rate.
    - `:wrat`: Surface water rate.
    - `:rate`: Total volumetric surface rate (upper limit).
    - `:rate_upper`: Total volumetric surface rate (upper limit).
    - `:rate_lower`: Total volumetric surface rate (lower limit).
    - `:resv`: Reservoir voidage.
- `val`: The value to which the limit is to be translated.

# Returns
- The translated limit value for the specified control parameter.
"""
function translate_limit(control::ProducerControl, name, val)
    # Note: Negative sign convention for production.
    # A lower absolute bound on a rate
    # |q| > |lim| -> q < lim if both sides are negative
    # means that we specify is_lower for upper limits and the other
    # way around for lower limits, when dealing with rates.
    is_lower = true
    if name == :bhp
        # Upper limit, pressure
        target_limit = BottomHolePressureTarget(val)
        # Pressures are positive, this is a true lower bound
        is_lower = true
    elseif name == :orat
        # Upper limit, surface oil rate
        target_limit = SurfaceOilRateTarget(val)
    elseif name == :lrat
        # Upper limit, surface liquid (water + oil) rate
        target_limit = SurfaceLiquidRateTarget(val)
    elseif name == :grat
        # Upper limit, surface gas rate
        target_limit = SurfaceGasRateTarget(val)
    elseif name == :wrat
        # Upper limit, surface water rate
        target_limit = SurfaceWaterRateTarget(val)
    elseif name == :rate || name == :rate_upper
        # Upper limit, total volumetric surface rate
        target_limit = TotalRateTarget(val)
    elseif name == :rate_lower
        # Lower limit, total volumetric surface rate. This is useful
        # disabling producers if they would otherwise start to inject.
        target_limit = TotalRateTarget(val)
        is_lower = false
    elseif name == :resv
        v, w = val
        target_limit = ReservoirVoidageTarget(v, w)
    else
        error("$name limit not supported for well acting as producer.")
    end
    return (target_limit, is_lower)
end

"""
    translate_limit(control::InjectorControl, name, val)

Translate the limit for a given `InjectorControl` object.

# Arguments
- `control::InjectorControl`: The control object for which the limit is being translated.
- `name`: The name of the limit to be translated.
    - `:bhp`: Bottom hole pressure.
    - `:rate`: Total volumetric surface rate (upper limit).
    - `:rate_upper`: Total volumetric surface rate (upper limit).
    - `:rate_lower`: Total volumetric surface rate (lower limit).
    - `:resv_rate`: Total volumetric reservoir rate.
- `val`: The value associated with the limit.

# Returns
- The translated limit value.
"""
function translate_limit(control::InjectorControl, name, val)
    is_lower = false
    if name == :bhp
        # Upper limit, pressure
        target_limit = BottomHolePressureTarget(val)
    elseif name == :rate || name == :rate_upper
        # Upper limit, total volumetric surface rate
        target_limit = TotalRateTarget(val)
    elseif name == :rate_lower
        # Lower limit, total volumetric surface rate
        target_limit = TotalRateTarget(val)
        is_lower = true
    elseif name == :resv_rate
        # Upper limit, total volumetric reservoir rate
        target_limit = TotalReservoirRateTarget(val)
    else
        error("$name limit not supported for well acting as injector.")
    end
    return (target_limit, is_lower)
end

function check_limit(current_control, target_limit, target, is_lower::Bool, q_t, source_model, well_state, rhoS, S)
    if typeof(target_limit) == typeof(target)
        # We are already operating at this target and there is no need to check.
        ok = true
        current_val = limit_val = NaN
    else
        current_val = value(well_target_value(q_t, current_control, target_limit, source_model, well_state, rhoS, S))
        limit_val = target_limit.value
        ϵ = 1e-6
        if is_lower
            # Limit is lower bound, check that we are above...
            ok = current_val >= (1 + ϵ)*limit_val
        else
            ok = current_val <= (1 - ϵ)*limit_val
        end
    end
    return (ok, current_val, limit_val)
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

"""
Well target contribution from well itself (disabled, zero value)
"""
function well_target(control, target::DisabledTarget, well_model, well_state, rhoS, S)
    return 0.0
end

"""
Well target contribution from well itself (bhp)
"""
function well_target(control, target::BottomHolePressureTarget, well_model, well_state, rhoS, S)
    return bottom_hole_pressure(well_state)
end

"""
Well target contribution from well itself (surface volume, injector)
"""
function well_target(control::InjectorControl, target::SurfaceVolumeTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    t_phases = lumped_phases(target)
    w_phases = get_phases(flow_system(well_model.system))
    t = 0.0
    for (ix, mix) in control.phases
        if w_phases[ix] in t_phases
            t += mix
        end
    end
    return t/control.mixture_density
end


"""
Well target contribution from well itself (reservoir volume, injector)
"""
function well_target(control::InjectorControl, target::TotalReservoirRateTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    w_phases = get_phases(flow_system(well_model.system))
    t = 0.0
    rho = well_state.PhaseMassDensities
    for (ix, mix) in control.phases
        t += mix/rho[ix, 1]
    end
    return t
end

"""
Well target contribution from well itself (surface volume, injector)
"""
function well_target(control::InjectorControl, target::TotalRateTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    if abs(target.value) == MIN_ACTIVE_WELL_RATE
        # Special meaning: Limits set at the CONST minimum value are really absolute mass limits.
        w = 1.0
    else
        w = 1.0/control.mixture_density
    end
    return w
end

"""
Well target contribution from well itself (surface volume, producer)
"""
function well_target(control::ProducerControl, target::SurfaceVolumeTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    phases = get_phases(flow_system(well_model.system))
    Tw = eltype(surface_volume_fractions)
    if abs(target.value) == MIN_ACTIVE_WELL_RATE
        # Special meaning: Limits set at the CONST minimum value are really absolute mass limits.
        w = 1.0
    else
        # Compute total density at surface conditions by weighting phase volumes at surf
        ρ_tot = zero(Tw)
        for (ρ, V) in zip(surface_densities, surface_volume_fractions)
            ρ_tot += ρ*V
        end
        # Divide by total density to get total volume at surface, then multiply that by surface volume fraction
        w = zero(Tw)
        if isa(target, TotalRateTarget)
            for i in eachindex(phases)
                @inbounds V = surface_volume_fractions[i]
                w += V
            end
        else
            lp = lumped_phases(target)
            for (i, ph) in enumerate(phases)
                if ph in lp
                    @inbounds V = surface_volume_fractions[i]
                    w += V
                end
            end
        end
        w = w/ρ_tot
    end
    return w
end

"""
Well target contribution from well itself (RESV, producer)
"""
function well_target(control::ProducerControl, target::ReservoirVoidageTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    Tw = eltype(surface_volume_fractions)
    ρ_tot = zero(Tw)
    for (ρ, V) in zip(surface_densities, surface_volume_fractions)
        ρ_tot += ρ*V
    end
    w = zero(Tw)
    for (i, S) in enumerate(surface_volume_fractions)
        w += S*target.weights[i]
    end
    return w/ρ_tot
end

function well_target(control::InjectorControl, target::ReinjectionTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    w = 1.0/control.mixture_density
    return w
end

function well_target_value(q_t, control, target, source_model, well_state, rhoS, S)
    v = well_target(control, target, source_model, well_state, rhoS, S)
    if rate_weighted(target)
        v *= value(q_t)
    end
    return v
end

