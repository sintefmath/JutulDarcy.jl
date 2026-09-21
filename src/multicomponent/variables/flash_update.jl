@jutul_secondary function update_flash!(flash_results, fr::FlashResults,
        model, Pressure, Temperature, OverallMoleFractions, ix)
    eos = model.system.equation_of_state
    p = Pressure
    T = Temperature
    z = OverallMoleFractions
    return flash_entity_loop!(flash_results, fr, eos, p, T, z, nothing, ix)
end

@jutul_secondary function update_flash!(flash_results, fr::FlashResults,
        model::LVCompositionalModel3Phase, Pressure, Temperature,
        OverallMoleFractions, ImmiscibleSaturation, ix)
    eos = model.system.equation_of_state
    flash_entity_loop!(flash_results, fr, eos, Pressure, Temperature,
        OverallMoleFractions, ImmiscibleSaturation, ix)
end

@inline immiscible_saturation(::Nothing, i) = 0.0
@inline immiscible_saturation(saturation, i) = @inbounds saturation[i]

@inline function flash_entity_loop!(flash_results, fr, eos, Pressure,
        Temperature, OverallMoleFractions, sw, ix)
    @inbounds for i in ix
        old = flash_results[i]
        flash_results[i] = immutable_flash_result(old, fr, eos,
            Pressure[i], Temperature[i], OverallMoleFractions,
            immiscible_saturation(sw, i), i)
    end
    return flash_results
end

@inline function estimated_K_from_previous_flash(f, V, z)
    # Estimate each component's vapor partition from the previous equilibrium,
    # then renormalize the two phases using the *current* overall composition.
    minimum_z = MultiComponentFlash.MINIMUM_COMPOSITION
    ratios = map(f.liquid.mole_fractions, f.vapor.mole_fractions) do x, y
        x = max(Float64(compositional_primal(x)), minimum_z)
        y = max(Float64(compositional_primal(y)), minimum_z)
        V*y/((1.0 - V)*x)
    end
    liquid_total = sum(map((ratio, z_i) -> z_i/(1.0 + ratio), ratios, z))
    vapor_total = sum(map((ratio, z_i) -> z_i*ratio/(1.0 + ratio), ratios, z))
    return ratios.*(liquid_total/vapor_total)
end

@inline function valid_reused_flash(V, K, converged)
    valid = converged && isfinite(V) && 1e-6 < V < 1.0 - 1e-6
    trivial = true
    @inbounds for K_i in K
        valid &= isfinite(K_i) && K_i > 0.0
        trivial &= abs(K_i - 1.0) < 1e-6
    end
    return valid && !trivial
end

@inline function numeric_flash(f,
        fr::FlashResults{M, use_stability_bypass, reuse_guess}, eos, cond
    ) where {M, use_stability_bypass, reuse_guess}
    config = MultiComponentFlash.StaticConfig()
    # This function may use approach by Rasmussen et al (2006) for reusing
    # previous flash results depending on options set
    if reuse_guess && f.state == MultiComponentFlash.two_phase_lv
        V0 = Float64(compositional_primal(f.V))
        if isfinite(V0) && 1e-6 < V0 < 1.0 - 1e-6
            K0 = estimated_K_from_previous_flash(f, V0, cond.z)
            # An interior vapor-fraction guess skips stability testing in
            # MultiComponentFlash. This is not guaranteed to be successful; it
            # must remain two-phase and nontrivial to be accepted.
            V, K, report = flash_2ph!(config, K0, eos, cond, V0;
                method = fr.method,
                maxiter = 20,
                extra_out = true,
                tolerance = fr.tolerance,
                z_min = nothing,
                stability_bypass = false,
                check = false,
                verbose = false)
            if valid_reused_flash(V, K, report.converged)
                return V, K, report.stability_result
            end
        end
    end

    # The quick attempt was disabled or failed. Start from Wilson K-values and
    # let the flash perform its stability test. The distance-based bypass is a
    # separate option. Two-phase history has a NaN distance, so it cannot bypass.
    K0 = initial_guess_K(eos, cond, config)
    if use_stability_bypass
        stability_storage = MultiComponentFlash.StaticStabilityStorage(
            f.flash_cond, f.critical_distance)
    else
        stability_storage = nothing
    end
    V, K, report = flash_2ph!(config, K0, eos, cond, NaN;
        method = fr.method,
        extra_out = true,
        tolerance = fr.tolerance,
        z_min = nothing,
        stability_storage = stability_storage,
        stability_bypass = use_stability_bypass,
        bypass_tolerance = fr.tolerance_bypass,
        check = false,
        verbose = false
    )
    return V, K, report.stability_result
end

@inline function pure_immiscible_flash(f, eos, cond)
    config = MultiComponentFlash.StaticConfig()
    K = initial_guess_K(eos, cond, config)
    report = MultiComponentFlash.StabilityReport(
        stable_liquid = true, stable_vapor = true)
    storage = MultiComponentFlash.StaticStabilityStorage(cond, NaN)
    stability = MultiComponentFlash.StaticStabilityResult(
        true, report, K, storage, false)
    return NaN, K, stability
end

@inline function immutable_flash_result(f,
        fr::FlashResults{M, use_stability_bypass},
        eos::GenericCubicEOS{E, R, N}, P, temperature,
        OverallMoleFractions, Sw, cell = 1) where {M, use_stability_bypass, E, R, N}
    z = cell_composition(Val(N), OverallMoleFractions, cell)
    z_numeric = numeric_composition(z)
    cond_numeric = (
        p = R(compositional_primal(P)),
        T = R(compositional_primal(temperature)),
        z = z_numeric)
    if is_pure_single_phase(compositional_primal(Sw))
        V_numeric, K_numeric, stability = pure_immiscible_flash(f, eos, cond_numeric)
    else
        V_numeric, K_numeric, stability = numeric_flash(f, fr, eos, cond_numeric)
    end

    Num = typeof(P + temperature + first(z))
    cond = (p = convert(Num, P), T = convert(Num, temperature), z = z)
    if isnan(V_numeric)
        is_vapor = single_phase_label(eos, cond_numeric) > 0.5
        state, V, x, y, Z_l, Z_v = single_phase_flash_result(
            eos, cond, cond_numeric, K_numeric, stability, is_vapor)
    else
        state, V, x, y, Z_l, Z_v = two_phase_flash_result(
            eos, cond, P, V_numeric, K_numeric, fr.tolerance)
    end

    K_out = K_numeric
    if use_stability_bypass
        cd = compositional_primal(stability.storage.critical_distance)
        critical_distance = Float64(cd)
        flash_cond = stability.storage.reference
    else
        critical_distance = NaN
        flash_cond = cond_numeric
    end
    return FlashedMixture2Phase(state, K_out, V, x, y, Z_l, Z_v,
        critical_distance, flash_cond, stability.report)
end

@inline function immutable_flash_result(f, fr,
        eos::KValuesEOS{E, R, N}, P, temperature,
        OverallMoleFractions, Sw, cell = 1) where {E, R, N}
    z = cell_composition(Val(N), OverallMoleFractions, cell)
    Num = typeof(P + temperature + first(z))
    z = SVector{N, Num}(z)
    cond = (p = convert(Num, P), T = convert(Num, temperature), z = z)
    K = SVector{N, Num}(initial_guess_K(eos, cond))
    V = MultiComponentFlash.solve_rachford_rice(K, z)
    if V <= zero(V)
        state = MultiComponentFlash.single_phase_l
        V = zero(V)
        x = y = z
    elseif V >= one(V)
        state = MultiComponentFlash.single_phase_v
        V = one(V)
        x = y = z
    else
        state = MultiComponentFlash.two_phase_lv
        x, y = phase_mole_fractions(z, K, V)
    end
    K_out = numeric_values(K)
    cond_numeric = (
            p = R(value(P)),
            T = R(value(temperature)),
            z = numeric_composition(z, R)
        )
    return FlashedMixture2Phase(state, K_out, V, x, y,
        one(Num), one(Num), NaN, cond_numeric, f.flash_stability)
end
