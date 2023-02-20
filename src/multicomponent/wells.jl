function flash_wellstream_at_surface(well_model, sys::S, well_state, rhoS) where S<:MultiComponentSystem
    total_masses = well_state.TotalMasses
    T = eltype(total_masses)
    nph = length(rhoS)
    rho = zeros(T, nph)
    volfrac = zeros(T, nph)
    eos = sys.equation_of_state
    nc = MultiComponentFlash.number_of_components(eos)

    if nph == 3
        a, l, v = phase_indices(sys)
        rhoWS = rhoS[a]
        # Convert to surface conditions
        rhoWW = well_state.PhaseMassDensities[a, 1]
        S_other = well_state.Saturations[a, 1]*rhoWW/rhoWS
        # Surface density given for aqueous phase
        rho[a] = rhoWS
        volfrac[a] = S_other
        n = nc + 1
    else
        l, v = phase_indices(sys)
        n = nc
        S_other = zero(T)
    end
    buf = InPlaceFlashBuffer(nc)

    sc = well_model.domain.grid.surface
    Pressure = sc.p
    Temperature = sc.T

    z = SVector{nc}(well_state.OverallMoleFractions[:, 1])
    m = SSIFlash()
    fstorage = flash_storage(eos, method = m, inc_jac = true, diff_externals = true, npartials = n, static_size = true)
    update_flash_buffer!(buf, eos, Pressure, Temperature, z)

    f = FlashedMixture2Phase(eos, T)
    x = f.liquid.mole_fractions
    y = f.vapor.mole_fractions
    forces = buf.forces

    result = update_flash_result(fstorage, m, eos, f.K, x, y, buf.z, forces, Pressure, Temperature, z)

    rho[l], rho[v] = mass_densities(eos, Pressure, Temperature, result)
    rem = one(T) - S_other
    S_l, S_v = phase_saturations(eos, Pressure, Temperature, result)

    volfrac[l] = rem*S_l
    volfrac[v] = rem*S_v
    rho = tuple(rho...)
    return (rho, volfrac)
end

Base.@propagate_inbounds function multisegment_well_perforation_flux!(out, sys::CompositionalSystem, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well

    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    X = state_res.LiquidMassFractions
    Y = state_res.VaporMassFractions

    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    X_w = state_well.LiquidMassFractions
    Y_w = state_well.VaporMassFractions

    nc = size(X, 1)
    nph = size(μ, 1)
    has_water = nph == 3
    phase_ix = phase_indices(sys)

    if has_water
        A, L, V = phase_ix
    else
        L, V = phase_ix
    end

    mob(ph) = kr[ph, rc]/μ[ph, rc]
    λ_t = zero(eltype(kr))
    for ph in 1:nph
        λ_t += mob(ph)
    end

    function phase_mass_flux(ph)
        dp = perforation_phase_potential_difference(conn, state_res, state_well, ph)
        # dp is pressure difference from reservoir to well. If it is negative,
        # we are injecting into the reservoir.
        if dp < 0
            mob_ph = λ_t*s_w[ph, wc]
            dens_ph = ρ_w[ph, wc]
        else
            mob_ph = mob(ph)
            dens_ph = ρ[ph, rc]
        end
        Q = mob_ph*dens_ph*dp
        return (Q, dp)
    end

    Q_l, dp_l = phase_mass_flux(L)
    Q_v, dp_v = phase_mass_flux(V)

    @inbounds for c in 1:nc
        if dp_l < 0
            # Injection
            X_upw = X_w[c, wc]
        else
            X_upw = X[c, rc]
        end
        if dp_v < 0
            Y_upw = Y_w[c, wc]
        else
            Y_upw = Y[c, rc]
        end
        out[c] = Q_l*X_upw + Q_v*Y_upw
    end
    if has_water
        out[nc+1], _ = phase_mass_flux(A)
    end
end

