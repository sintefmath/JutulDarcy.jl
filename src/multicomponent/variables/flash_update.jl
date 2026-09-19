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

@inline function numeric_flash(f, fr::FlashResults{M, StabilityBypass, false},
    eos, cond) where {M, StabilityBypass}
    return full_numeric_flash(f, fr, eos, cond)
end

@inline function numeric_flash(f,
        fr::FlashResults{M, StabilityBypass, true}, eos, cond) where {M, StabilityBypass}
    if f.state == MultiComponentFlash.two_phase_lv
        config = MultiComponentFlash.StaticConfig()
        V0 = Float64(compositional_primal(f.V))
        V, K, report = flash_2ph!(config, f.K, eos, cond, V0;
            method = SSIFlash(),
            maxiter = 20,
            extra_out = true,
            tolerance = fr.tolerance,
            z_min = nothing,
            stability_bypass = false,
            check = false,
            verbose = false)
        valid = report.converged && isfinite(V) &&
            1e-6 < V < 1.0 - 1e-6
        trivial = true
        @inbounds for K_i in K
            trivial &= abs(K_i - 1.0) < 1e-6
        end
        if valid && !trivial
            return V, K, report.stability_result
        end
    end
    return full_numeric_flash(f, fr, eos, cond)
end

@inline function full_numeric_flash(f,
        fr::FlashResults{M, StabilityBypass}, eos, cond) where {
        M, StabilityBypass}
    config = MultiComponentFlash.StaticConfig()
    K0 = initial_guess_K(eos, cond, config)
    storage = previous_stability_storage(f, Val(StabilityBypass))
    V, K, report = flash_2ph!(config, K0, eos, cond, NaN;
        method = SSIFlash(),
        extra_out = true,
        tolerance = fr.tolerance,
        z_min = nothing,
        stability_storage = storage,
        stability_bypass = StabilityBypass,
        bypass_tolerance = fr.tolerance_bypass,
        check = false,
        verbose = false)
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
        fr::FlashResults{M, StabilityBypass},
        eos::GenericCubicEOS{E, R, N}, P, temperature,
        OverallMoleFractions, Sw, cell = 1) where {
        M, StabilityBypass, E, R, N}
    z = cell_composition(Val(N), OverallMoleFractions, cell)
    z_numeric = numeric_composition(z)
    cond_numeric = (
        p = Float64(compositional_primal(P)),
        T = Float64(compositional_primal(temperature)),
        z = z_numeric)
    if is_pure_single_phase(compositional_primal(Sw))
        V_numeric, K_numeric, stability =
            pure_immiscible_flash(f, eos, cond_numeric)
    else
        V_numeric, K_numeric, stability =
            numeric_flash(f, fr, eos, cond_numeric)
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
    if StabilityBypass
        critical_distance = Float64(
            compositional_primal(stability.storage.critical_distance))
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
            p = Float64(value(P)),
            T = Float64(value(temperature)),
            z = numeric_composition(z)
        )
    return FlashedMixture2Phase(state, K_out, V, x, y,
        one(Num), one(Num), NaN, cond_numeric, f.flash_stability)
end
