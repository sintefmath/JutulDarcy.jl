function well_target_value(ctrl, target::ReservoirVoidageTarget, cond, well, model, state)
    q_w = cond.surface_aqueous_rate
    q_o = cond.surface_liquid_rate
    q_g = cond.surface_vapor_rate
    return compute_total_resv_rate(target.avg_state; qw = q_w, qo = q_o, qg = q_g)
end


function realize_control_for_reservoir(rstate, ctrl::ProducerControl{<:Union{HistoricalReservoirVoidageTarget, ReservoirVoidageTarget}}, model, dt)
    resv_t = ctrl.target
    state_avg = setup_average_resv_state(model, rstate)
    if resv_t isa HistoricalReservoirVoidageTarget
        qw = resv_t.water
        qo = resv_t.oil
        qg = resv_t.gas
        q_resv = compute_total_resv_rate(state_avg; qw = qw, qg = qg, qo = qo)
    else
        q_resv = resv_t.value
    end
    new_control = replace_target(ctrl, ReservoirVoidageTarget(q_resv, state_avg))
    return (new_control, true)
end

function setup_average_resv_state(model::StandardBlackOilModel, rstate; qw = 0.0, qo = 0.0, qg = 0.0)
    sys = model.system
    has_water = has_other_phase(sys)

    pv_t = 0.0
    p_avg = 0.0
    rs_avg = 0.0
    rv_avg = 0.0
    disgas = has_disgas(sys)
    vapoil = has_vapoil(sys)
    for c in eachindex(rstate.Pressure)
        p = value(rstate.Pressure[c])
        vol = value(rstate.FluidVolume[c])
        if has_water
            sw = value(rstate.ImmiscibleSaturation[c])
        else
            sw = 0.0
        end
        vol_hc = vol*(1.0 - sw)
        p_avg += p*vol_hc
        pv_t += vol_hc
        if disgas
            rs_avg += value(rstate.Rs[c])*vol_hc
        end
        if vapoil
            rv_avg += value(rstate.Rv[c])*vol_hc
        end
    end
    p_avg /= pv_t
    rs_avg /= pv_t
    rv_avg /= pv_t

    if has_water
        a, l, v = phase_indices(sys)
    else
        l, v = phase_indices(sys)
    end

    svar = Jutul.get_secondary_variables(model)
    b_var = svar[:ShrinkageFactors]
    reg = b_var.regions
    if has_water
        bW = shrinkage(b_var.pvt[a], reg, p_avg, 1)
    else
        bW = 1.0
    end
    if has_disgas(model.system)
        rs = min(rs_avg, sys.rs_max[1](p_avg))
        bO = rs -> shrinkage(b_var.pvt[l], reg, p_avg, rs, 1)
    else
        rs = 0.0
        bO = rs -> shrinkage(b_var.pvt[l], reg, p_avg, 1)
    end
    if has_vapoil(model.system)
        rv = min(rv_avg, sys.rv_max[1](p_avg))
        bG = rv -> shrinkage(b_var.pvt[v], reg, p_avg, rv, 1)
    else
        rv = 0.0
        bG = rv -> shrinkage(b_var.pvt[v], reg, p_avg, 1)
    end
    return (bW = bW, bO = bO, bG = bG, p = p_avg, rs = rs_avg, rv = rv_avg)
end


function compute_total_resv_rate(state_avg; qw = 0.0, qo = 0.0, qg = 0.0)
    if qw + qo + qg < 0.0
        sgn = -1.0
    else
        sgn = 1.0
    end
    # Surface rates
    qw = abs(qw)
    qo = abs(qo)
    qg = abs(qg)

    if qw <= 1e-20
        rs = 0.0
    else
        rs = min(qg/qo, state_avg.rs)
    end
    if qg <= 1e-20
        rv = 0.0
    else
        rv = min(qo/qg, state_avg.rv)
    end
    shrink = max(1.0 - rs*rv, 1e-20)
    # Water
    new_water_rate = qw/state_avg.bW
    # Oil
    bO = state_avg.bO(rs)
    new_oil_rate = (qo - rv*qg)/(bO*shrink)
    # Gas
    bG = state_avg.bG(rv)
    new_gas_rate = (qg - rs*qo)/(bG*shrink)
    resv_rate = sgn*(new_water_rate + new_oil_rate + new_gas_rate)
    return resv_rate
end
