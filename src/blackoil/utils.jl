function check_blackoil_pvt(model, p; rv = nothing, rs = nothing, reg = nothing, cell = 1)
    p = ForwardDiff.Dual(p, 1.0)
    ∂(x) = x.partials[1]
    ∂(x::Float64) = 0.0

    sys = model.system
    b = model.secondary_variables[:ShrinkageFactors]
    wat = has_other_phase(sys)
    pvto = b.pvt[1+wat]
    pvtg = b.pvt[2+wat]

    if has_disgas(sys)
        if isnothing(rs)
            rs = sys.rs_max(p)
        end
        B_o = 1/shrinkage(pvto, reg, p, rs, cell)
    else
        rs = 0.0
        B_o = 1/shrinkage(pvto, reg, p, cell)
    end

    if has_vapoil(sys)
        if isnothing(rv)
            rv = sys.rv_max(p)
        end
        B_g = 1/shrinkage(pvtg, reg, p, rv, cell)
    else
        rv = 0.0
        B_g = 1/shrinkage(pvtg, reg, p, cell)
    end
    return coats_pvt_tests(value(p), value(B_o), value(B_g), value(rv), value(rs), ∂(B_o), ∂(B_g), ∂(rv), ∂(rs))
end

function coats_pvt_tests(p, B_o, B_g, R_v, R_s, dB_o, dB_g, dRv, dRs)
    # Consistency checks for black oil PVT:
    # A Note on IMPES and Some IMPES-Based Simulation Models, K.H. Coats, SPE Journal, Vol. 5, No. 3, September 2000
    ok1 = (B_g - R_v*B_o)*dRs > dB_o*(1.0 - R_v*R_s)
    if !ok1
        @warn "First test failed: (B_g - R_v*B_o)*dRs < dB_o*(1.0 - R_v*R_s)" (B_g - R_v*B_o)*dRs dB_o*(1.0 - R_v*R_s) p B_o B_g R_v R_s dB_o dB_g dRv dRs
    end
    ok2 = (B_o - R_s*B_g)*dRv > dB_g*(1.0 - R_v*R_s)
    if !ok2
        @warn "Second test failed: (B_o - R_s*B_g)*dRv < dB_g*(1.0 - R_v*R_s)" (B_o - R_s*B_g)*dRv dB_g*(1.0 - R_v*R_s) p B_o B_g R_v R_s dB_o dB_g dRv dRs
    end
    ok3 = R_s * R_v < 1
    if !ok3
        @warn "Third test failed: R_s * R_v < 1" R_s * R_v p B_o B_g R_v R_s dB_o dB_g dRv dRs
    end
    ok4 = B_g - R_v * B_o > 0
    if !ok4
        @warn "Fourth test failed: B_g - R_v * B_o > 0" B_g - R_v * B_o p B_o B_g R_v R_s dB_o dB_g dRv dRs
    end
    ok5 = B_o - R_s * B_g > 0
    if !ok5
        @warn "Fifth test failed: B_o - R_s * B_g > 0" B_o - R_s * B_g p B_o B_g R_v R_s dB_o dB_g dRv dRs
    end
    return ok1 && ok2 && ok3 && ok4 && ok5
end
