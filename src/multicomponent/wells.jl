function flash_wellstream_at_surface(var, well_model, system::S, state, rhoS, cond = default_surface_cond()) where S<:MultiComponentSystem
    eos = system.equation_of_state
    nc = MultiComponentFlash.number_of_components(eos)
    z = SVector{nc}(state.OverallMoleFractions[:, 1])
    flash, _, surface_moles = get_separator_intermediate_storage(var, system, z, 2, 1)
    S_l, S_v, rho_l, rho_v = separator_flash(flash, system, surface_moles, eos, cond, z)
    return compositional_surface_densities(state, system, S_l, S_v, rho_l, rho_v)
end

function separator_flash(flash, system, surface_moles, eos, cond, z)
    result = separator_flash!(flash, eos, cond, z)
    rho_l, rho_v = mass_densities(eos, cond.p, cond.T, result)
    S_l, S_v = phase_saturations(eos, cond.p, cond.T, result)
    @assert S_l + S_v ≈ 1.0
    return (S_l, S_v, rho_l, rho_v)
end

function separator_flash(flash, system, surface_moles, eos::KValuesEOS, cond, z)
    if has_other_phase(system)
        _, rho_l, rho_v = reference_densities(system)
    else
        rho_l, rho_v = reference_densities(system)
    end
    # TODO Clean up, correct, derivatives?
    V = flash_2ph(eos, (p = cond.p, T = cond.T, z = z))
    S_l = 1.0 - V
    S_v = V
    return (S_l, S_v, rho_l, rho_v)
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

    λ_res = perforation_reservoir_mobilities(state_res, state_well, sys, rc, wc)
    mob(ph) = λ_res[ph]
    λ_t = sum(λ_res)

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
    return out
end


function separator_surface_flash!(var, model, system::MultiPhaseCompositionalSystemLV, state)
    eos = system.equation_of_state
    nc = MultiComponentFlash.number_of_components(eos)
    z0 = SVector{nc}(state.OverallMoleFractions[:, 1])
    rhoS, vol = separator_hydrocarbon_surface_flash!(var, system, z0, 1)
    return compositional_surface_densities(state, system, vol[1], vol[2], rhoS[1], rhoS[2])
end

function separator_hydrocarbon_surface_flash!(var, system::MultiPhaseCompositionalSystemLV, z0, well)
    # For each stage we need to keep track of total moles and composition.
    eos = system.equation_of_state
    F, moles, surface_moles = get_separator_intermediate_storage(var, system, z0, 2, well)
    moles[1] = z0
    separator_conditions = var.separator_conditions[well]
    separator_targets = var.separator_targets[well]
    n_stage = length(separator_conditions)
    for i in 1:n_stage
        cond = separator_conditions[i]
        streams, mole_fractions = flash_stream!(moles[i], F, eos, cond)
        targets = separator_targets[i]
        for ph in eachindex(targets)
            dest = targets[ph]
            q = streams[ph].*mole_fractions[ph]
            if dest == 0
                surface_moles[ph] += q
            else
                @assert dest > i "Destination was $dest for stage $i"
                moles[dest] += q
            end
        end
    end
    # Final stage: Flash the tank conditions for each component.
    cond = default_surface_cond()
    T = eltype(z0)
    nph = length(surface_moles)
    rhoS = @MVector zeros(T, nph)
    vol = @MVector zeros(T, nph)
    vol_total = zero(T)
    for ph in eachindex(surface_moles)
        sm = surface_moles[ph]
        tot_moles = sum(sm)
        if tot_moles ≈ 0
            vol[ph] = zero(T)
            rhoS[ph] = one(T)
        else
            z = sm./tot_moles
            result = separator_flash!(F, eos, cond, z)
            V = result.V
            L = one(V) - V
            rho_l, rho_v = mass_densities(eos, cond.p, cond.T, result)
            S_l, S_v = phase_saturations(eos, cond.p, cond.T, result)
            V_l = molar_volume(eos, cond.p, cond.T, result.liquid)
            V_v = molar_volume(eos, cond.p, cond.T, result.vapor)
            # This could in fact still be two-phase so we just sum up the
            # densities and volumes according to their fractions to get a volume
            # average and pretend it is single-phase. This way we preserve mass
            # in the effective density even if it is not strictly necessary.
            V_i = tot_moles*(L*V_l + V*V_v)
            vol[ph] = V_i
            vol_total += V_i
            # Volume weighted density
            rhoS[ph] = rho_l*S_l + rho_v*S_v
        end
    end
    # Normalize surface volumes to get fractions
    for ph in eachindex(vol)
        vol[ph] /= vol_total
    end
    return (rhoS, vol)
end

function separator_flash!(flash, eos, cond, z)
    settings = FlashResults{SSIFlash, false, false}(
        SSIFlash(), false, 1e-8, 10.0, false, false)
    return immutable_flash_result(
        flash, settings, eos, cond.p, cond.T, z, 0.0)
end

function flash_stream!(moles::SVector{N, T}, flash, eos, cond) where {N, T}
    total_moles = sum(moles)
    if total_moles ≈ 0
        q_l = q_v = zero(T)
        x = y = zero(flash.liquid.mole_fractions)
        result = flash
    else
        z = moles./total_moles
        result = separator_flash!(flash, eos, cond, z)
        V = result.V
        L = one(V) - V

        x = result.liquid.mole_fractions
        y = result.vapor.mole_fractions

        q_l = total_moles.*L
        q_v = total_moles.*V
    end

    return ((q_l, q_v), (x, y))
end

function get_separator_intermediate_storage(var, system::MultiPhaseCompositionalSystemLV, z::S, nph, well = 1) where S
    s = var.storage[well]
    z::Union{AbstractVector, Tuple, NamedTuple}
    T = eltype(z)
    if !haskey(s, T)
        eos = system.equation_of_state
        nc = MultiComponentFlash.number_of_components(eos)
        f = static_flashed_mixture(eos, T)
        n_stages = length(var.separator_conditions[well])
        # Array of vectors, one for each stage
        moles = zeros(S, n_stages)
        # Mutable surface moles
        surface_moles = zeros(S, nph)
        s[T] = (
            flash = f,
            moles = moles,
            surface_moles = surface_moles
        )
    end
    flash, moles, surface = s[T]
    @. moles = zero(S)
    @. surface = zero(S)
    return (flash, moles, surface)
end

function surface_density_and_volume_fractions(var, model, system::MultiPhaseCompositionalSystemLV, state, well)
    eos = system.equation_of_state
    n_hc = MultiComponentFlash.number_of_components(eos)
    mw = MultiComponentFlash.molar_masses(eos)
    rates = state.SurfaceComponentRates
    molar_rates = SVector{n_hc}(ntuple(c -> rates[c, well]/mw[c], Val(n_hc)))
    total_molar_rate = sum(molar_rates)
    if abs(value(total_molar_rate)) < MIN_ACTIVE_WELL_RATE
        z = map(_ -> one(total_molar_rate)/n_hc, molar_rates)
    else
        z = molar_rates/total_molar_rate
    end

    if isempty(var.separator_conditions[well])
        flash, _, surface_moles = get_separator_intermediate_storage(var, system, z, 2, well)
        S_l, S_v, rho_l, rho_v = separator_flash(
            flash, system, surface_moles, eos, default_surface_cond(), z)
    else
        rho, volume = separator_hydrocarbon_surface_flash!(var, system, z, well)
        rho_l, rho_v = rho
        S_l, S_v = volume
    end
    return compositional_surface_densities_from_rates(
        rates, well, system, S_l, S_v, rho_l, rho_v)
end

function compositional_surface_densities_from_rates(rates, well, system,
        S_l::S_T, S_v::S_T, rho_l::R_T, rho_v::R_T) where {S_T, R_T}
    nph = number_of_phases(system)
    T = promote_type(eltype(rates), S_T, R_T)
    if isbitstype(T)
        rhoS = @MVector zeros(T, nph)
        volume = @MVector zeros(T, nph)
    else
        rhoS = zeros(T, nph)
        volume = zeros(T, nph)
    end
    if has_other_phase(system)
        a, l, v = phase_indices(system)
        rhoWS = reference_densities(system)[a]
        n_hc = MultiComponentFlash.number_of_components(system.equation_of_state)
        q_hc = sum(c -> rates[c, well], 1:n_hc)
        q_water = rates[end, well]
        rho_hc = S_l*rho_l + S_v*rho_v
        volume_hc = q_hc/rho_hc
        volume_water = q_water/rhoWS
        total_volume = volume_hc + volume_water
        if abs(value(total_volume)) < MIN_ACTIVE_WELL_RATE
            water_fraction = zero(total_volume)
            hydrocarbon_fraction = one(total_volume)
        else
            water_fraction = volume_water/total_volume
            hydrocarbon_fraction = volume_hc/total_volume
        end
        rhoS[a] = rhoWS
        volume[a] = water_fraction
    else
        l, v = phase_indices(system)
        hydrocarbon_fraction = one(T)
    end
    rhoS[l] = rho_l
    rhoS[v] = rho_v
    volume[l] = S_l*hydrocarbon_fraction
    volume[v] = S_v*hydrocarbon_fraction
    if isbitstype(T)
        rhoS = SVector(rhoS)
        volume = SVector(volume)
    end
    return (rhoS, volume)
end

function compositional_surface_densities(state, system, S_l::S_T, S_v::S_T, rho_l::R_T, rho_v::R_T) where {S_T, R_T}
    nph = number_of_phases(system)
    T = promote_type(S_T, R_T)
    if isbitstype(T)
        rhoS = @MVector zeros(T, nph)
        volume = @MVector zeros(T, nph)
    else
        rhoS = zeros(T, nph)
        volume = zeros(T, nph)
    end
    if has_other_phase(system)
        a, l, v = phase_indices(system)
        rhoWS = reference_densities(system)[a]
        # Convert to surface conditions
        rhoWW = state.PhaseMassDensities[a, 1]
        S_other = state.Saturations[a, 1]*rhoWW/rhoWS
        # Surface density given for aqueous phase
        rhoS[a] = rhoWS
        volume[a] = S_other
        rem = one(T) - S_other
    else
        l, v = phase_indices(system)
        rem = one(T)
    end
    rhoS[l] = rho_l
    rhoS[v] = rho_v
    volume[l] = S_l*rem
    volume[v] = S_v*rem
    @assert sum(volume) ≈ 1.0 "Volume should sum to 1, was $(sum(volume))"
    if isbitstype(T)
        rhoS = SVector(rhoS)
        volume = SVector(volume)
    end
    return (rhoS, volume)
end
