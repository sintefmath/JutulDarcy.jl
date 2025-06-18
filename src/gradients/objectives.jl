"""
    compute_well_qoi(model::MultiModel, state, forces, well::Symbol, target::Union{WellTarget, Type})

Compute the quantity of interest (QoI) for a specified well in a reservoir simulation.

# Arguments
- `model::MultiModel`: The simulation model (from `setup_reservoir_model`).
- `state`: The current state of the simulation.
- `forces`: The forces applied in the simulation.
- `well::Symbol`: The symbol representing the well for which the QoI is computed.
- `target::Union{WellTarget, Type, Symbol}`: The target type (as type or symbol)
  or specific well target for the QoI computation.

Possible targets (as symbols):
- :gor (Gas-to-Oil Ratio)
- :wcut (Water Cut)
- :temperature (Well Temperature)
- :mass_rate (Total Surface Mass Rate)
- :bhp (Bottom Hole Pressure)
- :lrat (Surface Liquid Rate)
- :wrat (Surface Water Rate)
- :orat (Surface Oil Rate)
- :grat (Surface Gas Rate)

# Returns
- The computed QoI for the specified well.

"""
function compute_well_qoi(model::MultiModel, state, forces, well::Symbol, target::Union{WellTarget, Type, Symbol})
    well_model = model[well]
    rhoS = reference_densities(well_model.system)

    if haskey(model.models, :Facility)
        pos = get_well_position(model.models[:Facility].domain, well)
        fk = :Facility
    else
        pos = 1
        fk = Symbol("$(well)_ctrl")
    end
    ctrl = forces[fk].control[well]
    fstate = state[fk]
    wstate = state[well]

    if target isa Symbol
        if target == :gor
            grat = compute_well_qoi(model, state, forces, well, :grat)
            orat = compute_well_qoi(model, state, forces, well, :orat)
            return grat/max(orat, 1e-12)
        elseif target == :wcut
            wrat = compute_well_qoi(model, state, forces, well, :wrat)
            lrat = compute_well_qoi(model, state, forces, well, :lrat)
            return wrat/lrat
        elseif target == :temperature
            return wstate[:Temperature][1]
        elseif target == :mass_rate
            return fstate[:TotalSurfaceMassRate][pos]
        elseif target == :lrat
            target = SurfaceLiquidRateTarget
        elseif target == :wrat
            target = SurfaceWaterRateTarget
        elseif target == :orat
            target = SurfaceOilRateTarget
        elseif target == :grat
            target = SurfaceGasRateTarget
        elseif target == :bhp
            target = BottomHolePressureTarget
        else
            error("Unsupported QoI symbol $target.")
        end
    end

    translate_target_to_symbol
    if ctrl isa DisabledControl
        qoi = 0.0
    else
        if target isa Type
            if target<:SurfaceVolumeTarget
                if ctrl isa InjectorControl
                    tv = 1.0
                else
                    tv = -1.0
                end
            else
                tv = 100e5
            end
            target = target(tv)
        end
        if well_target_information(target).symbol != well_target_information(ctrl.target).symbol
            ctrl = replace_target(ctrl, target)
        end
        qoi = compute_well_qoi(well_model, wstate, fstate, well::Symbol, pos, rhoS, ctrl)
    end
    return qoi
end

function compute_well_qoi(well_model, well_state, fstate, well::Symbol, pos, rhoS, control)
    well_state = convert_to_immutable_storage(well_state)
    q_t = fstate[:TotalSurfaceMassRate][pos]
    target = control.target

    rhoS, S = surface_density_and_volume_fractions(well_state)
    v = well_target(control, target, well_model, well_state, rhoS, S)
    if rate_weighted(target)
        v *= q_t
    end
    return v
end


"""
    well_mismatch(qoi, wells, model_f, states_f, model_c, state_c, dt, step_info, forces; <keyword arguments>)

Compute well mismatch for a set of qoi's (well targets) and a set of well symbols.
"""
function well_mismatch(qoi, wells, model_f, states_f, model_c, state_c, dt, step_info, forces; weights = ones(length(qoi)), scale = 1.0, signs = nothing)
    if !(qoi isa AbstractArray)
        qoi = [qoi]
    end
    if !(wells isa AbstractArray)
        wells = [wells]
    end
    step_no = step_info[:step]
    obj = 0.0
    @assert length(weights) == length(qoi)
    for well in wells
        pos = get_well_position(model_c.models[:Facility].domain, well)

        well_f = model_f[well]
        well_c = model_c[well]
        rhoS = reference_densities(well_f.system)

        ctrl = forces[:Facility].control[well]
        if ctrl isa DisabledControl
            continue
        end

        state_f = states_f[step_no]

        for (i, q) in enumerate(qoi)
            ctrl = replace_target(ctrl, q)
            if !isnothing(signs)
                s = signs[i]
                if ctrl isa ProducerControl
                    sgn = -1
                else
                    sgn = 1
                end
                if s != sgn && s != 0
                    continue
                end
            end
            qoi_f = compute_well_qoi(well_f, state_f[well], state_f[:Facility], well, pos, rhoS, ctrl)
            qoi_c = compute_well_qoi(well_c, state_c[well], state_c[:Facility], well, pos, rhoS, ctrl)

            Δ = qoi_f - qoi_c
            obj += (weights[i]*Δ)^2
        end
    end
    return scale*dt*obj
end


"""
    setup_rate_optimization_objective(case, base_rate;
        max_rate_factor = 5,
        steps = :first,
        injectors = missing,
        producers = missing,
        verbose = true,
        limits_enabled = true,
        constraint = :total_sum_injected,
        kwarg...
    )

Setup the rate optimization objective and functions for a given case. The
objective is to maximize the NPV of the reservoir simulation by optimizing the
rates of the injectors. It is assumed that all injectors are initially set to
the same injection rate (`base_rate`). The largest value (if used in a
box-constrained setting) will be `max_rate_factor*base_rate`. If
`limits_enabled` is set to false, any well constraints will be ignored (other
than injectors switching to producers and vice versa).

Additional keyword arguments will be passed to the `npv_objective` function.
"""
function setup_rate_optimization_objective(case, base_rate;
        max_rate_factor = 5,
        steps = :first,
        injectors = missing,
        producers = missing,
        verbose = true,
        use_ministeps = true,
        limits_enabled = true,
        constraint = :total_sum_injected,
        sim_arg = NamedTuple(),
        kwarg...
    )
    steps in (:first, :each) || error("Invalid steps argument, must be :first or :each")
    function myprint(x)
        if verbose
            jutul_message("Rate optimization", x)
        end
    end
    max_rate = max_rate_factor*base_rate
    forces = case.forces
    eachstep = steps == :each
    nstep = length(case.dt)

    if forces isa Vector
        forces = [deepcopy(force) for force in forces]
    else
        forces = [deepcopy(forces) for i in 1:nstep]
    end

    case = JutulCase(case.model, case.dt, forces, state0 = case.state0, parameters = case.parameters)
    if eachstep
        myprint("steps=:$steps selected, optimizing controls for all $nstep steps separately.")
    else
        myprint("steps=:$steps selected, optimizing one set of controls for all steps.")
        for i in eachindex(case.forces)
            case.forces[i] = case.forces[1]
        end
    end
    if !limits_enabled
        for force in case.forces
            f_force = force[:Facility]
            for (ckey, ctrl) in f_force.control
                f_force.limits[ckey] = default_limits(ctrl)
            end
        end
    end
    wells = well_symbols(case.model)

    function find_wells(ctrl_t)
        out = Symbol[]
        for f in forces
            for w in wells
                ctrl = f[:Facility][:control][w]
                if ctrl isa ctrl_t
                    push!(out, w)
                end
            end
        end
        return unique(out)
    end
    if ismissing(injectors)
        injectors = find_wells(InjectorControl)
    end
    if ismissing(producers)
        producers = find_wells(ProducerControl)
    end
    ninj = length(injectors)

    length(injectors) > 0 || error("No injectors found")
    length(producers) > 0 || error("No producers found")
    is_both = intersect(injectors, producers)
    length(is_both) == 0 || error("Wells were both producers and injectors in forces? $is_both")
    myprint("$ninj injectors and $(length(producers)) producers selected.")
    cache = Dict()

    function f!(x; grad = true)
        nstep = length(case.dt)
        if eachstep
            nstep_unique = nstep
        else
            nstep_unique = 1
        end
        x = reshape(x, ninj, nstep_unique)
        for stepno in eachindex(case.forces)
            for (inj_no, inj) in enumerate(injectors)
                if eachstep
                    x_i = x[inj_no, stepno]
                else
                    x_i = x[inj_no, 1]
                end
                new_rate = max(x_i*max_rate, MIN_INITIAL_WELL_RATE)
                f_forces = case.forces[stepno][:Facility]
                ctrl = f_forces.control[inj]
                lims = f_forces.limits[inj]
                new_target = TotalRateTarget(new_rate)
                f_forces.control[inj] = replace_target(ctrl, new_target)
                f_forces.limits[inj] = merge(lims, as_limit(new_target))
            end
        end
        simulated = simulate_reservoir(case; output_substates = use_ministeps, info_level = -1, sim_arg...)
        r = simulated.result
        dt_mini = report_timesteps(r.reports, ministeps = use_ministeps)
        function npv_obj(model, state, dt, step_info, forces)
            return npv_objective(model, state, dt, step_info, forces;
                injectors = injectors,
                producers = producers,
                timesteps = dt_mini,
                kwarg...
            )
        end
        obj = Jutul.evaluate_objective(npv_obj, case.model, r.states, case.dt, case.forces)
        targets = Jutul.force_targets(case.model)
        targets[:Facility][:control] = :control
        targets[:Facility][:limits] = nothing
        if grad
            if !haskey(cache, :storage)
                cache[:storage] = Jutul.setup_adjoint_forces_storage(case.model, r.states, forces, case.dt, npv_obj;
                    state0 = case.state0,
                    targets = targets,
                    parameters = case.parameters,
                    eachstep = eachstep,
                    di_sparse = true
                )
            end
            dforces, t_to_f, grad_adj = Jutul.solve_adjoint_forces!(cache[:storage], case.model, r.states, r.reports, npv_obj, forces,
                state0 = case.state0,
                parameters = case.parameters
            )
            df = zeros(ninj, nstep_unique)
            for stepno in 1:nstep_unique
                for (inj_no, inj) in enumerate(injectors)
                    ctrl = dforces[stepno][:Facility].control[inj]
                    do_dq = ctrl.target.value
                    df[inj_no, stepno] = do_dq*max_rate
                end
            end
            out = (obj, vec(df))
        else
            out = obj
        end
        return out
    end

    if constraint == :total_sum_injected
        if eachstep
            # TODO: This requires updates in lbfgs
            # A = spzeros(nstep, ninj*nstep)
            A = zeros(nstep, ninj*nstep)
            for i in 1:nstep
                offset = (i-1)*ninj
                A[i, (offset+1):(offset+ninj)] .= 1
            end
            b = fill(ninj*base_rate/max_rate, nstep)
        else
            A = ones(1, length(injectors))
            b = [ninj*base_rate/max_rate]
        end
        lin_eq = (A=A, b=b)
    else
        @assert constraint == :none || isnothing(constraint) "Constraint must be :total_sum_injected or :none"
        lin_eq = NamedTuple()
    end
    if eachstep
        x0 = fill(base_rate/max_rate, ninj*nstep)
    else
        x0 = fill(base_rate/max_rate, ninj)
    end
    # Get just objective
    function F!(x)
        return f!(x, grad = false)
    end
    # Get objective and update gradient in-place
    function F_and_dF!(dFdx, x)
        obj, grad = f!(x, grad = true)
        dFdx .= grad
        return obj
    end
    # Get just the gradient
    function dF!(dFdx, x)
        obj, grad = f!(x, grad = true)
        dFdx .= grad
        return dFdx
    end

    return (
        x0 = x0,
        lin_eq = lin_eq,
        obj = f!,
        F_and_dF! = F_and_dF!,
        dF! = dF!,
        F! = F!,
        case = case
    )
end


"""
    npv_objective(model, state, dt, step_info, forces;
        timesteps,
        injectors,
        producers,
        oil_price = 60.0,
        gas_price = 10.0,
        water_price = -3.0,
        water_cost = 5.0,
        oil_cost = oil_price,
        gas_cost = gas_price,
        liquid_unit = si_unit(:stb),
        gas_unit = si_unit(:kilo)*si_unit(:feet)^3,
        discount_rate = 0.025,
        discount_unit = si_unit(:year),
        scale = 1.0
    )

Evaluate the contribution to net-present-value for a given step, with the given
costs and prices for oil, gas, and water and discount rate. Costs are assumed to
be the cost of injecting a given fluid and prices are the revenue from producing
the corresponding fluids. Prices and costs can be set to negative (e.g. to
account for water being a cost when produced).

Setting `maximize` to false will make the objective function negative-valued,
with larger negative values corresponding to better results. This is useful for
some optimizers.
"""
function npv_objective(model, state, dt, step_info, forces;
        timesteps,
        injectors,
        producers,
        maximize = true,
        oil_price = 60.0,
        gas_price = 10.0,
        water_price = -3.0,
        water_cost = 5.0,
        oil_cost = oil_price,
        gas_cost = gas_price,
        liquid_unit = si_unit(:stb),
        gas_unit = si_unit(:kilo)*si_unit(:feet)^3,
        discount_rate = 0.025,
        discount_unit = si_unit(:year),
        scale = 1.0
    )
    step_no = step_info[:step]
    phases = get_phases(reservoir_model(model).system)
    has_wat = AqueousPhase() in phases
    has_gas = VaporPhase() in phases
    has_oil = LiquidPhase() in phases

    obj = 0.0
    for w in producers
        if has_oil
            orat = -compute_well_qoi(model, state, forces, w, SurfaceOilRateTarget)/liquid_unit
        else
            orat = 0.0
        end
        if has_gas
            grat = -compute_well_qoi(model, state, forces, w, SurfaceGasRateTarget)/gas_unit
        else
            grat = 0.0
        end
        if has_wat
            wrat = -compute_well_qoi(model, state, forces, w, SurfaceWaterRateTarget)/liquid_unit
        else
            wrat = 0.0
        end
        obj += oil_price*orat + gas_price*grat + water_price*wrat
    end
    for w in injectors
        if has_oil
            orat = compute_well_qoi(model, state, forces, w, SurfaceOilRateTarget)/liquid_unit
        else
            orat = 0.0
        end
        if has_gas
            grat = compute_well_qoi(model, state, forces, w, SurfaceGasRateTarget)/gas_unit
        else
            grat = 0.0
        end
        if has_wat
            wrat = compute_well_qoi(model, state, forces, w, SurfaceWaterRateTarget)/liquid_unit
        else
            wrat = 0.0
        end
        obj -= (oil_cost*orat + gas_cost*grat + water_cost*wrat)
    end
    time = step_info[:time] + dt
    if maximize
        sgn = 1
    else
        sgn = -1
    end

    return sgn*dt*obj*((1.0+discount_rate)^(-time/discount_unit))/scale
end
