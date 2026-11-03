# Rasmussen et al. (2006), Eqs. (5)-(6); Belkadi et al. (2011), Sec. 2.1.
# Keep component distribution fractions, not the old total vapor fraction.
@inline function estimated_phase_amounts_from_previous_flash(f, V, z)
    F = eltype(z)
    V = F(V)
    liquid = map(f.liquid.mole_fractions, f.vapor.mole_fractions, z) do x, y, zi
        l = (one(V) - V) * F(compositional_primal(x))
        v = V * F(compositional_primal(y))
        zi * (l / (l + v))
    end
    vapor = map(f.liquid.mole_fractions, f.vapor.mole_fractions, z) do x, y, zi
        l = (one(V) - V) * F(compositional_primal(x))
        v = V * F(compositional_primal(y))
        zi * (v / (l + v))
    end
    # Compute both amounts directly: z_i - v_i loses the trace liquid amount
    # when a component is almost entirely in the vapor phase (or vice versa).
    return liquid, vapor
end

@inline function reused_flash_is_stable(eos, cond, K, V, tolerance)
    # Phase-dependent models can have additional stationary splits. Equality
    # of fugacities alone cannot distinguish a metastable split from equilibrium.
    MultiComponentFlash.forces_per_phase(eos) || return true
    F = eltype(K)
    x, _ = phase_mole_fractions(cond.z, K, V)
    forces = MultiComponentFlash.static_force_coefficients(eos, cond, F)
    phase_cond = (p = cond.p, T = cond.T, z = x, phase = Val(:liquid))
    reference = MultiComponentFlash.static_fugacities(eos, phase_cond, forces, F)
    K0 = initial_guess_K(eos, cond, MultiComponentFlash.StaticConfig())
    return reused_flash_trial_is_stable(
        eos, cond, reference, K0, forces,
        Val(true), tolerance
    ) && reused_flash_trial_is_stable(
        eos, cond,
        reference, K0, forces, Val(false), tolerance
    )
end

@inline function reused_flash_trial_is_stable(
        eos, cond, reference, K0, forces,
        inside_is_vapor::Val, tolerance
    )
    stable, _, _, _, trial = MultiComponentFlash.static_michelsen_test(
        reference, cond.z, K0, eos, cond, forces, inside_is_vapor
    )
    trial_cond = (
        p = cond.p, T = cond.T, z = trial,
        phase = MultiComponentFlash.stability_phase(inside_is_vapor),
    )
    fugacity = MultiComponentFlash.static_fugacities(
        eos, trial_cond, forces, eltype(K0)
    )
    tpd = sum(trial .* log.(fugacity ./ reference))
    return stable && isfinite(tpd) && tpd >= -tolerance
end

@inline function valid_reused_flash(V, K, converged, z, tolerance)
    valid = converged && isfinite(V) && 1.0e-6 < V < 1.0 - 1.0e-6
    valid &= all(k -> isfinite(k) && k > 0, K)
    valid || return false
    x, y = phase_mole_fractions(z, K, V)
    # A trace component can have a huge K even when the phases are numerically
    # indistinguishable. Such a split cannot support a reliable linearization.
    return maximum(abs, x - y) > tolerance
end
