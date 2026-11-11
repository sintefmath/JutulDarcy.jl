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

@inline function valid_reused_flash(eos, V, K, converged, z, tolerance)
    valid = converged && isfinite(V) && 1.0e-6 < V < 1.0 - 1.0e-6
    valid &= all(k -> isfinite(k) && k > 0, K)
    valid || return false
    # A reversed SW split can also be stationary with the wrong correlations.
    if MultiComponentFlash.forces_per_phase(eos)
        @inbounds K[gas_water_component_index(eos)] < 1 || return false
    end
    x, y = phase_mole_fractions(z, K, V)
    # A trace component can have a huge K even when the phases are numerically
    # indistinguishable. Such a split cannot support a reliable linearization.
    return maximum(abs, x - y) > tolerance
end
