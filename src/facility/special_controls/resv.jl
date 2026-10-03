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

    disgas = has_disgas(sys)
    vapoil = has_vapoil(sys)
    Tv = eltype(rstate.FluidVolume)
    TP = eltype(rstate.Pressure)
    Ts = eltype(rstate.ImmiscibleSaturation)
    TRs = eltype(rstate.Rs)
    TPv = eltype(rstate.Rv)

    T = promote_type(Tv, TP, Ts, TRs, TPv)
    if T <: Jutul.AdjointsDI.SparseConnectivityTracer.Dual
        # Hack for SparseConnectivityTracer.GradientTracer
        return (
            bW = 1.0,
            bO = _ -> 1.0,
            bG = _ -> 1.0,
            p = 100e5,
            rs = 0.0,
            rv = 0.0
        )
    end
    if has_water
        hc_pv_fn = Base.Broadcast.broadcasted(
            (vol, sw) -> value(vol)*(1.0 - value(sw)),
            rstate.FluidVolume, rstate.ImmiscibleSaturation)
        pv_t = sum(hc_pv_fn)
    else
        pv_t = sum(value, rstate.FluidVolume)
    end
    function weighted_average(x)
        if has_water
            wfn = Base.Broadcast.broadcasted(
                (v, vol, sw) -> value(v)*value(vol)*(1.0 - value(sw)),
                x, rstate.FluidVolume, rstate.ImmiscibleSaturation)
        else
            wfn = Base.Broadcast.broadcasted(
                (v, vol) -> value(v)*value(vol),
                x, rstate.FluidVolume)
        end
        return sum(wfn)/pv_t
    end
    p_avg = weighted_average(rstate.Pressure)
    if disgas
        rs_avg = weighted_average(rstate.Rs)
    else
        rs_avg = 0.0
    end
    if vapoil
        rv_avg = weighted_average(rstate.Rv)
    else
        rv_avg = 0.0
    end

    svar = Jutul.get_secondary_variables(model)
    b_var = svar[:ShrinkageFactors]
    reg = b_var.regions
    if has_water
        a, l, v = phase_indices(sys)
        bW = shrinkage(b_var.pvt[a], reg, p_avg, 1)
    else
        l, v = phase_indices(sys)
        bW = 1.0
    end
    if disgas
        rs_max = sys.rs_max[1](p_avg)
        bO = rs -> shrinkage(b_var.pvt[l], reg, p_avg, min(rs, rs_max), 1)
    else
        rs_max = 0.0
        bO = rs -> shrinkage(b_var.pvt[l], reg, p_avg, 1)
    end
    if vapoil
        rv_max = sys.rv_max[1](p_avg)
        bG = rv -> shrinkage(b_var.pvt[v], reg, p_avg, min(rv, rv_max), 1)
    else
        rv_max = 0.0
        bG = rv -> shrinkage(b_var.pvt[v], reg, p_avg, 1)
    end
    return (bW = bW, bO = bO, bG = bG, p = p_avg, rs = rs_avg, rv = rv_avg)
end

function compute_total_resv_rate(state_avg; qw = 0.0, qo = 0.0, qg = 0.0, is_obs::Bool = false)
    if qw + qo + qg < 0.0
        sgn = -1.0
    else
        sgn = 1.0
    end
    # Surface rates
    qw = abs(qw)
    qo = abs(qo)
    qg = abs(qg)

    rs = clamp(qg/(qo + 1e-12), 0.0, state_avg.rs)
    rv = clamp(qo/(qg + 1e-12), 0.0, state_avg.rv)
    shrink = max(1.0 - rs*rv, 1e-20)
    bW = state_avg.bW
    bO = state_avg.bO(rs)
    bG = state_avg.bG(rv)

    # Water
    new_water_rate = qw/bW
    # Oil
    new_oil_rate = (qo - rv*qg)/(bO*shrink)
    # Gas
    new_gas_rate = (qg - rs*qo)/(bG*shrink)
    resv_rate = new_water_rate + new_oil_rate + new_gas_rate

    return sgn*resv_rate
end
