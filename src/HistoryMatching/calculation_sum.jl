function (local_obj::SumHistoryMatchObjective)(model, state, dt, step_info, forces)
    local_obj.evaluation_count[] += 1
    hm = local_obj.match
    ctrls = forces[:Facility].control
    fmodel = model[:Facility]
    fstate = state[:Facility]

    val = get_step_mismatch_contributions(hm, fmodel, fstate, ctrls, step_info)
    return val*hm.total_scale
end


function get_step_mismatch_contributions(hm, fmodel, fstate, ctrls, step_info)
    val = 0.0
    eval_match(x, sgn, target) = squared_mismatch_for_step(hm.logger, x, fmodel, fstate, ctrls, sgn, Val(target), step_info)
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

function squared_mismatch_for_step(logger, matches::Vector{WellMatch}, fmodel, fstate, ctrls, sgn, ::Val{target}, step_info) where target
    out = 0.0
    for wm in matches
        wname = wm.name
        val, obs, w_for_step = mismatch_for_step(fmodel, fstate, ctrls[wname], sgn, wm, target, step_info, missing)
        dt = step_info[:dt]
        squared_delta = w_for_step*(val - obs)^2
        added_value = dt*squared_delta
        if !ismissing(logger.data)
            start = stop = step_info[:substep_global]
            if !haskey(logger.data, wname)
                logger.data[wname] = []
            end
            tsym = JutulDarcy.translate_target_to_symbol(target)
            push!(logger.data[wname], (start = start, stop = stop, target = tsym, value = added_value, is_cumulative = false))
        end
        out += added_value
    end
    return out
end
