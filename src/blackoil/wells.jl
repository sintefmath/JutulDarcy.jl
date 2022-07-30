function well_perforation_flux!(out, sys::BlackOilSystem, state_res, state_well, rhoS, WI, rc, wc)
    p_res = state_res.Pressure
    p_well = state_well.Pressure


    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    rs = state_res.Rs
    b = state_res.ShrinkageFactors
    # s = state_res.Saturations

    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    rs_w = state_well.Rs
    b_w = state_well.ShrinkageFactors
    # TODO: Fix the pressure drop
    ρgdz = 0
    a, l, v = 1, 2, 3
    dp = WI*(p_well[wc] - p_res[rc] + ρgdz)
    λ_a = kr[a, rc]/μ[a, rc]
    λ_l = kr[l, rc]/μ[l, rc]
    λ_v = kr[v, rc]/μ[v, rc]
    rhoOS = rhoS[l]
    rhoGS = rhoS[v]

    if dp >= 0
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
    out[a] = -q_a
    out[l] = -q_l
    out[v] = -q_v
end

function apply_well_reservoir_sources!(sys::BlackOilSystem, res_q, well_q, state_res, state_well, param_res, param_well, perforations, sgn)
    p_res = state_res.Pressure
    p_well = state_well.Pressure

    val = x -> local_ad(x, nothing)

    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    rs = state_res.Rs
    b = state_res.ShrinkageFactors
    s = state_res.Saturations

    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    rs_w = state_well.Rs
    b_w = state_well.ShrinkageFactors

    rhoS = param_res[:reference_densities]

    perforation_sources_blackoil!(well_q, perforations, val(p_res),         p_well,  val(kr), val(s), val(μ), val(ρ), val(b), val(rs),    ρ_w,      b_w,      s_w,      rs_w,  sgn, rhoS)
    perforation_sources_blackoil!(res_q,  perforations,     p_res,      val(p_well),     kr,      s,     μ,      ρ,      b,      rs, val(ρ_w), val(b_w), val(s_w), val(rs_w), sgn, rhoS)
end

function perforation_sources_blackoil!(target, perf, p_res, p_well, kr, s, μ, ρ, b, rs, ρ_w, b_w, s_w, rs_w, sgn, rhoS)
    # (self -> local cells, reservoir -> reservoir cells, WI -> connection factor)
    nc = size(ρ, 1)
    nph = size(μ, 1)
    a, l, v = 1, 2, 3

    rhoOS = rhoS[l]
    rhoGS = rhoS[v]
    @inbounds for i in eachindex(perf.self)
        si, ri, wi, gdz = unpack_perf(perf, i)
        if gdz != 0
            dens_w = @views mix_by_saturations(s_w[:, si], ρ_w[:, si])
            dens_r = @views mix_by_saturations(s[:, ri], ρ[:, ri])
            ρ_mix = 0.5*(dens_w + dens_r)
            ρgdz = gdz*ρ_mix
        else
            ρgdz = 0
        end
        @inbounds dp = wi*(p_well[si] - p_res[ri] + ρgdz)
        λ_a = kr[a, ri]/μ[a, ri]
        λ_l = kr[l, ri]/μ[l, ri]
        λ_v = kr[v, ri]/μ[v, ri]
        if dp >= 0
            # Injection
            λ_t = λ_a + λ_l + λ_v
            Q = sgn*λ_t*dp

            bO = b_w[l, si]
            bG = b_w[v, si]
            rs_i = rs_w[si]

            sO = s_w[l, si]
            sG = s_w[v, si]

            target[a, i] = s_w[a, si]*ρ_w[a, si]*Q
            target[l, i] = sO*rhoOS*bO*Q
            target[v, i] = rhoGS*(sG*bG + rs_i*sO*bO)*Q
        else
            # Production
            Q = sgn*dp
            target[a, i] = Q*ρ[a, ri]*λ_a

            bO = b[l, ri]
            bG = b[v, ri]

            α_l = bO*λ_l
            α_v = bG*λ_v

            q_o = Q*bO*λ_l*rhoOS
            q_g = Q*(bO*λ_l*rs[ri] + bG*λ_v)*rhoGS

            target[l, i] = q_o
            target[v, i] = q_g
        end
    end
end


function flash_wellstream_at_surface(well_model::SimulationModel{D, S}, well_state, rhoS) where {D, S<:BlackOilSystem}
    vol = well_state.TotalMasses[:, 1]./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end
