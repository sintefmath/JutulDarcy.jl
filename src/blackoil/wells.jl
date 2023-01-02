Base.@propagate_inbounds function well_perforation_flux!(out, sys::StandardBlackOilSystem, state_res, state_well, rhoS, dp, rc, wc)
    λ_a, λ_l, λ_v = res_mobility(state_res, sys, rc)
    a, l, v, rhoGS, rhoOS = well_pvt_bo(sys)
    b, b_w, ρ, ρ_w, s_w = well_volumes_bo(state_res, state_well)
    rs = state_res.Rs
    rv = state_res.Rv
    rs_w = state_well.Rs
    rv_w = state_well.Rv
    if dp < 0.0
        # Injection
        λ_t = λ_a + λ_l + λ_v
        Q = λ_t*dp
        bO = b_w[l, wc]
        bG = b_w[v, wc]
        rs_i = rs_w[wc]
        rv_i = rv_w[wc]
        sO = s_w[l, wc]
        sG = s_w[v, wc]
        q_a = s_w[a, wc]*ρ_w[a, wc]*Q
        q_l = rhoOS*(sO*bO + rv_i*sG*bG)*Q
        q_v = rhoGS*(sG*bG + rs_i*sO*bO)*Q
    else
        # Production
        bO = b[l, rc]
        bG = b[v, rc]
        q_a = dp*ρ[a, rc]*λ_a
        q_l = dp*(bO*λ_l + bG*λ_v*rv[rc])*rhoOS
        q_v = dp*(bO*λ_l*rs[rc] + bG*λ_v)*rhoGS
    end
    out[a] = q_a
    out[l] = q_l
    out[v] = q_v
end

Base.@propagate_inbounds function well_perforation_flux!(out, sys::VapoilBlackOilSystem, state_res, state_well, rhoS, dp, rc, wc)
    λ_a, λ_l, λ_v = res_mobility(state_res, sys, rc)
    a, l, v, rhoGS, rhoOS = well_pvt_bo(sys)
    b, b_w, ρ, ρ_w, s_w = well_volumes_bo(state_res, state_well)
    rv = state_res.Rv
    rv_w = state_well.Rv
    if dp < 0.0
        # Injection
        λ_t = λ_a + λ_l + λ_v
        Q = λ_t*dp
        bO = b_w[l, wc]
        bG = b_w[v, wc]
        rv_i = rv_w[wc]
        sO = s_w[l, wc]
        sG = s_w[v, wc]
        q_a = s_w[a, wc]*ρ_w[a, wc]*Q
        q_l = rhoOS*(sO*bO + rv_i*sG*bG)*Q
        q_v = rhoGS*sG*bG*Q
    else
        # Production
        bO = b[l, rc]
        bG = b[v, rc]
        q_a = dp*ρ[a, rc]*λ_a
        q_l = dp*(bO*λ_l + bG*λ_v*rv[rc])*rhoOS
        q_v = dp*bG*λ_v*rhoGS
    end
    out[a] = q_a
    out[l] = q_l
    out[v] = q_v
end

Base.@propagate_inbounds function well_perforation_flux!(out, sys::DisgasBlackOilSystem, state_res, state_well, rhoS, dp, rc, wc)
    λ_a, λ_l, λ_v = res_mobility(state_res, sys, rc)
    a, l, v, rhoGS, rhoOS = well_pvt_bo(sys)
    b, b_w, ρ, ρ_w, s_w = well_volumes_bo(state_res, state_well)
    rs = state_res.Rs
    rs_w = state_well.Rs
    if dp < 0.0
        # Injection
        λ_t = λ_a + λ_l + λ_v
        Q = λ_t*dp
        bO = b_w[l, wc]
        bG = b_w[v, wc]
        rs_i = rs_w[wc]
        sO = s_w[l, wc]
        sG = s_w[v, wc]
        q_a = s_w[a, wc]*ρ_w[a, wc]*Q
        q_l = rhoOS*sO*bO*Q
        q_v = rhoGS*(sG*bG + rs_i*sO*bO)*Q
    else
        # Production
        bO = b[l, rc]
        bG = b[v, rc]
        q_a = dp*ρ[a, rc]*λ_a
        q_l = dp*bO*λ_l*rhoOS
        q_v = dp*(bO*λ_l*rs[rc] + bG*λ_v)*rhoGS
    end
    out[a] = q_a
    out[l] = q_l
    out[v] = q_v
end

function flash_wellstream_at_surface(well_model, system::S, well_state, rhoS) where S<:BlackOilSystem
    vol = well_state.TotalMasses[:, 1]./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end

function well_pvt_bo(sys)
    a, l, v = phase_indices(sys)
    rhoS = reference_densities(sys)
    rhoOS = rhoS[l]
    rhoGS = rhoS[v]
    return (a, l, v, rhoGS, rhoOS)
end

function res_mobility(state_res, sys, rc)
    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    a, l, v = phase_indices(sys)
    λ_a = kr[a, rc]/μ[a, rc]
    λ_l = kr[l, rc]/μ[l, rc]
    λ_v = kr[v, rc]/μ[v, rc]
    return (λ_a, λ_l, λ_v)
end

function well_volumes_bo(state_res, state_well)
    ρ = state_res.PhaseMassDensities
    b = state_res.ShrinkageFactors
    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    b_w = state_well.ShrinkageFactors
    return (b, b_w, ρ, ρ_w, s_w)
end
