Base.@propagate_inbounds function well_perforation_flux!(out, sys::BlackOilSystem, state_res, state_well, rhoS, dp, rc, wc)
    # Res properties
    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    rs = state_res.Rs
    b = state_res.ShrinkageFactors
    # Well properties
    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    rs_w = state_well.Rs
    b_w = state_well.ShrinkageFactors
    a, l, v = 1, 2, 3
    λ_a = kr[a, rc]/μ[a, rc]
    λ_l = kr[l, rc]/μ[l, rc]
    λ_v = kr[v, rc]/μ[v, rc]
    rhoOS = rhoS[l]
    rhoGS = rhoS[v]

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
        q_l = sO*rhoOS*bO*Q
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

function flash_wellstream_at_surface(well_model::SimulationModel{D, S}, well_state, rhoS) where {D, S<:BlackOilSystem}
    vol = well_state.TotalMasses[:, 1]./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end
