export TotalSurfaceMassRate, WellGroup, DisabledControl
export HistoryMode, PredictionMode, Wells

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

function Jutul.count_entities(wg::WellGroup, ::Wells)
    return length(wg.well_symbols)
end

function Jutul.get_domain_intersection(u::Cells, target_d::DiscretizedDomain{W}, source_d::WellControllerDomain,
                                           target_symbol, source_symbol) where {W<:WellGrid}
    # From controller to top well cell
    pos = get_well_position(source_d, target_symbol)
    if isnothing(pos)
        t = nothing
        s = nothing
    else
        t = [1]
        s = [pos]
    end
    (target = t, source = s, target_entity = u, source_entity = Wells())
end

function Jutul.get_domain_intersection(u::Wells, target_d::WellControllerDomain, source_d::DiscretizedDomain{W},
                                           target_symbol, source_symbol) where {W<:WellGrid}
    # From top cell in well to control equation
    pos = get_well_position(target_d, source_symbol)
    if isnothing(pos)
        t = nothing
        s = nothing
    else
        t = [pos]
        s = [1]
    end
    (target = t, source = s, target_entity = u, source_entity = Cells())
end

function get_well_position(d, symbol)
    match = findall(d.well_symbols .== symbol)
    if length(match) == 0
        return nothing
    else
        return only(match)
    end
end

function Jutul.associated_entity(::TotalSurfaceMassRate) Wells() end

function Jutul.update_primary_variable!(state, massrate::TotalSurfaceMassRate, state_symbol, model, dx)
    v = state[state_symbol]
    symbols = model.domain.well_symbols
    cfg = state.WellGroupConfiguration
    # Injectors can only have strictly positive injection rates,
    # producers can only have strictly negative and disabled controls give zero rate.
    function do_update(v, dx, ctrl)
        return Jutul.update_value(v, dx)
    end
    function do_update(v, dx, ctrl::InjectorControl)
        return Jutul.update_value(v, dx, nothing, nothing, MIN_ACTIVE_WELL_RATE, nothing)
    end
    function do_update(v, dx, ctrl::ProducerControl)
        return Jutul.update_value(v, dx, nothing, nothing, nothing, -MIN_ACTIVE_WELL_RATE)
    end
    function do_update(v, dx, ctrl::DisabledControl)
        # Set value to zero since we know it is correct.
        return Jutul.update_value(v, -value(v))
    end
    @inbounds for i in eachindex(v)
        s = symbols[i]
        v[i] = do_update(v[i], dx[i], operating_control(cfg, s))
    end
end


## Well controls

"""
Impact from well group in facility on conservation equation inside well
"""
function Jutul.update_cross_term!(ct::InjectiveCrossTerm, eq::ConservationLaw, well_storage, facility_storage,
                            target_model::SimulationModel{D, S}, source_model::SimulationModel{WG},
                            well_symbol, source, dt) where
                            {D<:DiscretizedDomain{W} where W<:WellGrid,
                            S<:MultiPhaseSystem,
                            WG<:WellGroup}
    fstate = facility_storage.state
    wstate = well_storage.state
    # Stuff from facility
    mswell = source_model.domain
    pos = get_well_position(mswell, well_symbol)
    cfg = fstate.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    qT = fstate.TotalSurfaceMassRate[pos]

    if isa(ctrl, InjectorControl)
        if value(qT) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(target_model.system)
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(qT) > 0
            @warn "Producer $well_symbol is injecting?"
        end
        top_node = 1
        masses = wstate.TotalMasses[:, top_node]
        mass = sum(masses)
        mix = masses./mass
    end
    update_topnode_sources!(ct.crossterm_source, ct.crossterm_target, qT, mix)
end

function update_topnode_sources!(cts, ctt, qT, mix)
    @inbounds for i in eachindex(mix)
        m = -mix[i]
        cts[i] = value(m)*qT
        ctt[i] = m*value(qT)
    end
end

"""
Cross term from well on control equation for well
"""
function Jutul.update_cross_term!(ct::InjectiveCrossTerm, eq::ControlEquationWell,
                            target_storage, source_storage,
                            target_model::SimulationModel{WG},
                            source_model::SimulationModel{D},
                            target, well_symbol, dt) where {D<:DiscretizedDomain{W} where W<:WellGrid, WG<:WellGroup}
    fstate = target_storage.state
    cfg = fstate.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    # q_t = facility_surface_mass_rate_for_well(target_model, well_symbol, fstate)
    well_state = source_storage.state
    param = source_storage.parameters
    rhoS = param[:reference_densities]

    update_facility_control_crossterm!(ct.crossterm_source, ct.crossterm_target, well_state, rhoS, target_model, source_model, ctrl, well_symbol, fstate)
end

function update_facility_control_crossterm!(s_buf, t_buf, well_state, rhoS, target_model, source_model, ctrl, well_symbol, fstate)
    target = ctrl.target
    q_t = facility_surface_mass_rate_for_well(target_model, well_symbol, fstate)
    if isa(target, DisabledTarget)
        # Early return - no cross term needed.
        t_∂w = value(q_t)
        t_∂f = q_t
        t_num = 0.0
    else
        cfg = fstate.WellGroupConfiguration
        limits = current_limits(cfg, well_symbol)

        is_injecting = value(q_t) >= 0
        has_limits = !isnothing(limits)
        is_bhp = isa(target, BottomHolePressureTarget)

        need_rates = isa(ctrl, ProducerControl) && (!is_bhp || has_limits)
        if need_rates
            rhoS, S = flash_wellstream_at_surface(source_model, well_state, rhoS)
        else
            S = nothing
        end
        if has_limits
            target = apply_well_limit!(cfg, target, source_model, well_state, well_symbol, rhoS, S, value(q_t), limits)
        end
        # Compute target value with AD relative to well.
        t = well_target(ctrl, target, source_model, well_state, rhoS, S)
        t_∂w = t
        t_∂f = value(t)

        if rate_weighted(target)
            t_∂w *= value(q_t)
            t_∂f *= q_t
        end
        t_num = target.value
    end
    scale = target_scaling(target)
    s_buf[1] = (t_∂w - t_num)/scale
    t_buf[1] = (t_∂f - t_num)/scale
end

rate_weighted(t) = true
rate_weighted(::BottomHolePressureTarget) = false
rate_weighted(::DisabledTarget) = false

target_scaling(::Any) = 1.0
target_scaling(::BottomHolePressureTarget) = 1e5

function facility_surface_mass_rate_for_well(model::SimulationModel, wsym, fstate)
    pos = get_well_position(model.domain, wsym)
    return fstate.TotalSurfaceMassRate[pos]
end

bottom_hole_pressure(ws) = ws.Pressure[1]

function surface_target_phases(target::SurfaceVolumeTarget, phases)
    return findall(in(lumped_phases(target)), phases)
end

surface_target_phases(target::TotalRateTarget, phases) = eachindex(phases)

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
    return 1.0/control.mixture_density
end

"""
Well target contribution from well itself (surface volume, producer)
"""
function well_target(control::ProducerControl, target::SurfaceVolumeTarget, well_model, well_state, surface_densities, surface_volume_fractions)
    phases = get_phases(well_model.system)
    positions = surface_target_phases(target, phases)

    @assert length(positions) > 0
    Tw = eltype(surface_volume_fractions)
    # Compute total density at surface conditions by weighting phase volumes at surf
    ρ_tot = zero(Tw)
    for (ρ, V) in zip(surface_densities, surface_volume_fractions)
        ρ_tot += ρ*V
    end
    # Divide by total density to get total volume at surface, then multiply that by surface volume fraction
    w = zero(Tw)
    for pos in positions
        V = surface_volume_fractions[pos]
        w += V
    end
    w = w/ρ_tot
    return w
end

function well_target_value(q_t, control, target, arg...)
    v = well_target(control, target, arg...)
    if rate_weighted(target)
        v *= value(q_t)
    end
    return v
end


function Jutul.associated_entity(::ControlEquationWell) Wells() end

function Jutul.update_equation!(eq::ControlEquationWell, storage, model, dt)
    state = storage.state
    ctrl = state.WellGroupConfiguration.operating_controls
    wells = model.domain.well_symbols
    surf_rate = state.TotalSurfaceMassRate
    entries = eq.equation.entries
    control_equations!(entries, wells, ctrl, surf_rate)
end

function control_equations!(entries, wells, ctrl, surf_rate)
    @inbounds for (i, key) in enumerate(wells)
        entries[i] = 0.0
    end
end

function Jutul.align_to_jacobian!(eq::ControlEquationWell, jac, model, u::Cells; kwarg...)
    # Need to align to cells, faces is automatically done since it is on the diagonal bands
    cache = eq.equation_top_cell
    layout = matrix_layout(model.context)
    control_equation_top_cell_alignment!(cache, jac, layout; kwarg...)
end

function control_equation_top_cell_alignment!(cache, jac, layout; equation_offset = 0, variable_offset = 0)
    nu, ne, np = ad_dims(cache)
    cellix = 1
    for e in 1:ne
        for d = 1:np
            pos = find_jac_position(jac, 1 + equation_offset, cellix + variable_offset, e, d, nu, nu, ne, np, layout)
            set_jacobian_pos!(cache, 1, e, d, pos)
        end
    end
end

# Selection of primary variables
function select_primary_variables_domain!(S, domain::WellGroup, system, formulation)
    S[:TotalSurfaceMassRate] = TotalSurfaceMassRate()
end

function select_equations_domain!(eqs, domain::WellGroup, system, arg...)
    # eqs[:potential_balance] = (PotentialDropBalanceWell, 1)
    eqs[:control_equation] = (ControlEquationWell, 1)
end

function setup_forces(model::SimulationModel{D}; control = nothing, limits = nothing, set_default_limits = true) where {D <: WellGroup}
    # error() # Fix me. Set up defaults for all wells, including rate limits if not provided.
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

function Jutul.initialize_extra_state_fields_domain!(state, model, domain::WellGroup)
    # Insert structure that holds well control (limits etc) that is then updated before each step
    state[:WellGroupConfiguration] = WellGroupConfiguration(domain.well_symbols)
end

function Jutul.update_before_step_domain!(storage, model::SimulationModel, domain::WellGroup, dt, forces)
    # Set control to whatever is on the forces
    cfg = storage.state.WellGroupConfiguration
    q_t = storage.state.TotalSurfaceMassRate
    op_ctrls = cfg.operating_controls
    req_ctrls = cfg.requested_controls
    for key in keys(forces.control)
        # If the requested control in forces differ from the one we are presently using, we need to switch.
        # Otherwise, stay the course.
        newctrl = forces.control[key]
        oldctrl = req_ctrls[key]
        if newctrl != oldctrl
            # We have a new control. Any previous control change is invalid.
            # Set both operating and requested control to the new one.
            @debug "Well $key switching from $oldctrl to $newctrl"
            req_ctrls[key] = newctrl
            op_ctrls[key] = newctrl
        end
        pos = get_well_position(model.domain, key)
        q_t[pos] = valid_surface_rate_for_control(q_t[pos], newctrl)
    end
    for key in keys(forces.limits)
        cfg.limits[key] = forces.limits[key]
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
        target, changed, current_val, limit_val, lim_type = check_active_limits(ctrl, target, current_lims, wmodel, wstate, well, density_s, volume_fraction_s, total_mass_rate)
        if changed
            old = cfg.operating_controls[well].target
            next = replace_target(ctrl, target)
            cfg.operating_controls[well] = next
            @debug "$well: Switching control from $old to $target due to $(typeof(target)) limit:\nComputed value $current_val exceeds $lim_type limit $limit_val.\nNew control: $next"
        end
    end
    return target
end

function check_active_limits(control, target, limits, wmodel, wstate, well::Symbol, density_s, volume_fraction_s, total_mass_rate)
    changed = false
    cval = tval = NaN
    is_lower = false
    for (name, val) in pairs(limits)
        if isfinite(val)
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
    else
        error("$name limit not supported for well acting as producer.")
    end
    return (target_limit, is_lower)
end

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
    else
        error("$name limit not supported for well acting as producer.")
    end
    return (target_limit, is_lower)
end

function check_limit(current_control, target_limit, target, is_lower::Bool, q_t, arg...)
    if typeof(target_limit) == typeof(target)
        # We are already operating at this target and there is no need to check.
        ok = true
        current_val = limit_val = NaN
    else
        current_val = value(well_target_value(q_t, current_control, target_limit, arg...))
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


function convergence_criterion(model, storage, eq::ControlEquationWell, r; dt = 1)
    wells = model.domain.well_symbols
    cfg = storage.state.WellGroupConfiguration
    e = abs.(vec(r))
    names = map(w -> name_equation(w, cfg), wells)
    R = Dict("Abs" => (errors = e, names = names))
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
