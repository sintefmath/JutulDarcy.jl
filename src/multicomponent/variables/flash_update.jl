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
