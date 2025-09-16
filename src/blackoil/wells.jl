Base.@propagate_inbounds function multisegment_well_perforation_flux!(out, sys::StandardBlackOilSystem, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well
    if has_other_phase(sys)
        a, l, v = phase_indices(sys)
        dp_a, dp_l, dp_v = res_dp(conn, state_res, state_well, sys)
        λ_a, λ_l, λ_v = perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc)
        λ_t = λ_a + λ_l + λ_v
        a, l, v, rhoGS, rhoOS = well_pvt_bo(sys)
        b, b_w, ρ, ρ_w, s_w = well_volumes_bo(state_res, state_well)
        # Water component flux
        if dp_a < 0.0
            # Injection
            Q_a = s_w[a, wc]*ρ_w[a, wc]*λ_t*dp_a
        else
            # Production
            Q_a = ρ[a, rc]*λ_a*dp_a
        end
        out[a] = Q_a
    else
        l, v = phase_indices(sys)
        dp_l, dp_v = res_dp(conn, state_res, state_well, sys)
        λ_l, λ_v = perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc)
        λ_t = λ_l + λ_v
        l, v, rhoGS, rhoOS = well_pvt_bo_2ph(sys)
        b, b_w, ρ, ρ_w, s_w = well_volumes_bo(state_res, state_well)
    end

    q_l = q_v = zero(dp_l)
    # Oil component flux
    if dp_l < 0.0
        # Injection
        bO = b_w[l, wc]
        sO = s_w[l, wc]
        q = λ_t*dp_l*sO*bO
        q_l += q
        if has_disgas(sys)
            q_v += state_well.Rs[wc]*q
        end
    else
        # Production
        bO = b[l, rc]
        q = dp_l*bO*λ_l
        q_l += q
        if has_disgas(sys)
            q_v += state_res.Rs[rc]*q
        end
    end
    # Gas component flux
    if dp_v < 0.0
        # Injection
        bG = b_w[v, wc]
        sG = s_w[v, wc]
        q = λ_t*dp_v*sG*bG
        q_v += q
        if has_vapoil(sys)
            q_l += state_well.Rv[wc]*q
        end
    else
        # Production
        bG = b[v, rc]
        q = dp_v*bG*λ_v
        q_v += q
        if has_vapoil(sys)
            q_l += state_res.Rv[rc]*q
        end
    end
    out[l] = q_l*rhoOS
    out[v] = q_v*rhoGS
    return out
end

Base.@propagate_inbounds function simple_well_perforation_flux!(out, sys::StandardBlackOilSystem, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well
    Q_a = Q_v = Q_l = 0

    Q_in = 0
    ρ = state_res.PhaseMassDensities
    if has_other_phase(sys)
        a, l, v = phase_indices(sys)
        dp_a, dp_l, dp_v = res_dp(conn, state_res, state_well, sys)
        λ_a, λ_l, λ_v = perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc)
        λ_t = λ_a + λ_l + λ_v

        a, l, v, rhoGS, rhoOS = well_pvt_bo(sys)
        # Water component flux
        if dp_a < 0.0
            # Injection
            Q_in += λ_a*ρ[a, rc]*dp_a
        else
            # Production
            Q_a += ρ[a, rc]*λ_a*dp_a
        end
    else
        l, v = phase_indices(sys)
        dp_l, dp_v = res_dp(conn, state_res, state_well, sys)
        λ_l, λ_v = perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc)
        λ_t = λ_l + λ_v
        l, v, rhoGS, rhoOS = well_pvt_bo_2ph(sys)
    end
    ρ = state_res.PhaseMassDensities
    b = state_res.ShrinkageFactors

    # Oil component flux
    if dp_l < 0.0
        # Injection
        Q_in += λ_l*ρ[l, rc]*dp_l
    else
        # Production
        bO = b[l, rc]
        q = dp_l*bO*λ_l
        Q_l += rhoOS*q
        if has_disgas(sys)
            Q_v += rhoGS*state_res.Rs[rc]*q
        end
    end
    # Gas component flux
    if dp_v < 0.0
        # Injection
        Q_in += λ_v*ρ[v, rc]*dp_v
    else
        # Production
        bG = b[v, rc]
        q = dp_v*bG*λ_v
        Q_v += rhoGS*q
        if has_vapoil(sys)
            Q_l += rhoOS*state_res.Rv[rc]*q
        end
    end

    if Q_in < 0.0
        X = state_well.MassFractions
        if has_other_phase(sys)
            Q_a += X[a]*Q_in
        end
        Q_l += X[l]*Q_in
        Q_v += X[v]*Q_in
    end
    if has_other_phase(sys)
        out[a] = Q_a
    end
    out[l] = Q_l
    out[v] = Q_v
end

function flash_wellstream_at_surface(var, well_model, system::S, well_state, rhoS, cond = default_surface_cond()) where S<:BlackOilSystem
    if haskey(well_state, :MassFractions)
        X = well_state.MassFractions
    else
        X = well_state.TotalMasses[:, 1]
    end
    vol = X./rhoS
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

function well_pvt_bo_2ph(sys)
    l, v = phase_indices(sys)
    rhoS = reference_densities(sys)
    rhoOS = rhoS[l]
    rhoGS = rhoS[v]
    return (l, v, rhoGS, rhoOS)
end


function perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc)
    λ = state_res.PhaseMobilities
    phases = phase_indices(sys)
    if haskey(state_res, :FullyMixedPolymerViscosityMultiplier)
        water = phases[1]
        function F(ph)
            m = λ[ph, rc]
            if ph == water
                m0 = m
                old = state_res.FullyMixedPolymerViscosityMultiplier[1, rc]
                pure = state_well.FullyMixedPolymerViscosityMultiplier[1, wc]
                m *= old/pure
            end
            return m
        end

    else
        F = ph -> λ[ph, rc]
    end
    mob = map(F, phases)
    return mob
end


function res_dp(conn, state_res, state_well, sys)
    return map(
        x -> perforation_phase_potential_difference(conn, state_res, state_well, x),
        phase_indices(sys))
end


function well_volumes_bo(state_res, state_well)
    ρ = state_res.PhaseMassDensities
    b = state_res.ShrinkageFactors
    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    b_w = state_well.ShrinkageFactors
    return (b, b_w, ρ, ρ_w, s_w)
end
