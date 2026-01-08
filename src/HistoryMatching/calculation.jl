function evaluate_global_match(hm::HistoryMatch, result::ReservoirSimResult)
    c = hm.case
    model = c.model
    states = result.states
    timesteps = c.dt
    all_forces = c.forces
    return Jutul.evaluate_objective(hm, model, result.result.states, timesteps, all_forces)
end

function get_well_value(wm::WellMatch, ctrl, target, fmodel, fstate)
    wname = wm.name
    val = JutulDarcy.compute_well_qoi(fmodel, fstate, wname, ctrl, target)
    return val
end

function get_well_observation(wm::WellMatch, time)
    return wm.data(time)
end

function get_well_match(wm::WellMatch, ctrl, target, fmodel, fstate, report_step, time)
    qoi = get_well_value(wm, ctrl, target, fmodel, fstate)
    obs = get_well_observation(wm, time)
    if (target isa BottomHolePressureTarget && obs < 1e-10) || isnan(obs)
        # Missing data
        obs = 0.0
        qoi = 0.0
    end
    w = effective_weight(wm, report_step)
    return (w*qoi, w*obs)
end

function (global_obj::GlobalHistoryMatchObjective)(model, state0, states, step_infos, forces, input_data)
    global_obj.evaluation_count[] += 1
    hm = global_obj.match
    obj = 0.0
    # Find start and stop index
    # println("Evaluate global history match objective #$(global_obj.evaluation_count[])")
    period_contribution(start_idx, stop_idx, weights) = get_period_contribution(hm, model, states, step_infos, forces, start_idx, stop_idx, weights)
    if ismissing(hm.periods)
        for i in eachindex(step_infos)
            obj += period_contribution(i, i, missing)
        end
    else
        for period in hm.periods
            period::MatchPeriod
            start_report_idx = first(period.step_idx)
            stop_report_idx = last(period.step_idx)
            start_idx = findfirst(x -> x[:step] == start_report_idx, step_infos)
            stop_idx = findlast(x -> x[:step] == stop_report_idx, step_infos)
            obj += period_contribution(start_idx, stop_idx, period.weights)
        end
    end
    obj += get_cumulative_contribution(hm, model, states, step_infos, forces)
    return obj*hm.total_scale
end

function (local_obj::SumHistoryMatchObjective)(model, state, dt, step_info, forces)
    local_obj.evaluation_count[] += 1
    hm = local_obj.match
    val = get_period_contribution(hm, model, (state,), (step_info,), (forces,), 1, 1, missing)
    return val*hm.total_scale
end

function get_cumulative_contribution(hm::HistoryMatch, model, states, step_infos, forces)
    start_idx = 1
    stop_idx = length(step_infos)
    eval_match(x, sgn, target) = get_period_contribution_wells(hm.logger, x, sgn, target, hm.wellpos, model, states, step_infos, forces, start_idx, stop_idx, missing, is_cumulative = true)
    prod_sgn = -1.0
    val = 0.0
    # Producers only...
    val += eval_match(hm.producer_cumulative_gas, prod_sgn, SurfaceGasRateTarget(-1.0))
    val += eval_match(hm.producer_cumulative_oil, prod_sgn, SurfaceOilRateTarget(-1.0))
    val += eval_match(hm.producer_cumulative_liquid, prod_sgn, SurfaceLiquidRateTarget(-1.0))
    val += eval_match(hm.producer_cumulative_water, prod_sgn, SurfaceWaterRateTarget(-1.0))
    return val
end

function get_period_contribution(hm::HistoryMatch, model, states, step_infos, forces, start_idx::Int, stop_idx::Int, weights)
    val = 0.0
    eval_match(x, sgn, target) = get_period_contribution_wells(hm.logger, x, sgn, target, hm.wellpos, model, states, step_infos, forces, start_idx, stop_idx, weights)
    inj_sgn = 1.0
    prod_sgn = -1.0
    val = 0.0
    # Injectors
    bhp_inj = eval_match(hm.injector_bhp, 1.0, BottomHolePressureTarget(1.0))
    rate_inj = eval_match(hm.injector_rate, inj_sgn, TotalRateTarget(1.0))
    orat_inj = eval_match(hm.injector_orat, inj_sgn, SurfaceOilRateTarget(1.0))
    wrat_inj = eval_match(hm.injector_wrat, inj_sgn, SurfaceWaterRateTarget(1.0))
    grat_inj = eval_match(hm.injector_grat, inj_sgn, SurfaceGasRateTarget(1.0))
    # Producers
    bhp_prod = eval_match(hm.producer_bhp, 1.0, BottomHolePressureTarget(1.0))
    rate_prod = eval_match(hm.producer_rate, prod_sgn, TotalRateTarget(-1.0))
    grat_prod = eval_match(hm.producer_grat, prod_sgn, SurfaceGasRateTarget(-1.0))
    orat_prod = eval_match(hm.producer_orat, prod_sgn, SurfaceOilRateTarget(-1.0))
    lrat_prod = eval_match(hm.producer_lrat, prod_sgn, SurfaceLiquidRateTarget(-1.0))
    wrat_prod = eval_match(hm.producer_wrat, prod_sgn, SurfaceWaterRateTarget(-1.0))

    if false
        println("Well match contributions: ")
        println(" BHP inj: $bhp_inj, rate inj: $rate_inj, orat inj: $orat_inj, wrat inj: $wrat_inj, grat inj: $grat_inj | ")
        println(" BHP prod: $bhp_prod, rate prod: $rate_prod, grat prod: $grat_prod, orat prod: $orat_prod, lrat prod: $lrat_prod, wrat prod: $wrat_prod ")
    end
    val = bhp_inj + rate_inj + orat_inj + wrat_inj + grat_inj +
        bhp_prod + rate_prod + grat_prod + orat_prod + lrat_prod + wrat_prod
    return val
end

function get_period_contribution_wells(logger::HistoryMatchLogger, wms::Vector{WellMatch}, sgn, target::JutulDarcy.WellTarget, wellpos, model, states, step_infos, forces, start_idx::Int, stop_idx::Int, weights; kwarg...)
    val = 0.0
    if length(wms) == 0
        return val
    end
    for wm in wms
        val += get_period_contribution_well(logger, wm, wellpos, sgn, target, model, states, step_infos, forces, start_idx, stop_idx, weights; kwarg...)
    end
    return val
end

function get_period_contribution_well(logger::HistoryMatchLogger, wm::WellMatch, wellpos, sgn, target::JutulDarcy.WellTarget, model, states, step_infos, forces, start_idx::Int, stop_idx::Int, weights;
        is_cumulative::Bool = false
    )
    wname = wm.name
    fmodel = model[:Facility]
    w_total = 0.0
    if is_cumulative
        @assert JutulDarcy.rate_weighted(target)
        obj = 0.0
        calculated_cumulative = 0.0

        for step_index in start_idx:stop_idx
            step_info = step_infos[step_index]
            control_step = step_info[:step]::Int
            w_step_from_period = ismissing(weights) ? 1.0 : weights[control_step]
            if w_step_from_period > 0.0
                dt = step_info[:dt]::Float64
                time = step_info[:time]::Float64
                # Unpack stuff
                state = states[step_index]
                force = forces[step_index]
                fstate = state[:Facility]
                ctrl = force[:Facility].control[wname]
                val = get_well_value(wm, ctrl, target, fmodel, fstate)
                # Note: Observations are cumulative, we need to calculate the volumes ourselves.
                calculated_cumulative += sgn*val*dt
                observed_cumulative = get_well_observation(wm, time)
                w_well = effective_weight(wm, control_step)
                # Accumulate objective
                obj += ((calculated_cumulative - observed_cumulative)*w_step_from_period*w_well)^2
            end
        end
        out = obj
    else
        calculated = 0.0
        observed = 0.0
        for step_index in start_idx:stop_idx
            step_info = step_infos[step_index]
            control_step = step_info[:step]::Int
            w_step_from_period = ismissing(weights) ? 1.0 : weights[control_step]
            if w_step_from_period > 0.0
                dt = step_info[:dt]::Float64
                time = step_info[:time]::Float64
                if JutulDarcy.rate_weighted(target)
                    w_dt = w_step_from_period*dt
                else
                    w_dt = w_step_from_period
                end
                # Unpack stuff
                state = states[step_index]
                force = forces[step_index]
                fstate = state[:Facility]
                ctrl = force[:Facility].control[wname]
                val, obs = get_well_match(wm, ctrl, target, fmodel, fstate, control_step, time)
                calculated += sgn*val*w_dt
                observed += obs*w_dt
                w_total += w_step_from_period
            end
        end
        out = ((calculated - observed)/w_total)^2
    end
    if !ismissing(logger.data)
        start = step_infos[1][:substep_global]
        stop = step_infos[end][:substep_global]
        logger.data[(well = wname, start = start, stop = stop, target = target)] = out
    end
    return out
end

function effective_weight(wm::WellMatch, step::Int)
    if wm.weight isa AbstractVector
        w = wm.weight[step]
    else
        w = wm.weight
    end
    return w*wm.scale
end

