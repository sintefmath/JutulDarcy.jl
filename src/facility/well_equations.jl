# local_discretization(e::PotentialDropBalanceWell, i) = e.flow_discretization(i, Faces())

Jutul.discretization(e::PotentialDropBalanceWell) = e.flow_discretization

function Jutul.update_equation_in_entity!(eq_buf, i, state, state0, eq::PotentialDropBalanceWell, model, dt, ldisc = local_discretization(eq, i))
    (; face, left, right, gdz) = ldisc
    μ = state.PhaseViscosities
    V = state.TotalMassFlux[face]
    densities = state.PhaseMassDensities
    s = state.Saturations
    p = state.Pressure

    rho_l, mu_l = saturation_mixed(s, densities, μ, left)
    rho_r, mu_r = saturation_mixed(s, densities, μ, right)

    Δθ = Jutul.two_point_potential_drop(p[left], p[right], gdz, rho_l, rho_r)
    if Δθ > 0
        μ_mix = mu_l
    else
        μ_mix = mu_r
    end
    rho = 0.5*(rho_l + rho_r)
    μ_mix = 0.5*(mu_l + mu_r)

    seg_model = model.domain.grid.segment_models[face]
    Δp = segment_pressure_drop(seg_model, V, rho, μ_mix)

    eq_buf[] = Δθ + Δp
end

function saturation_mixed(saturations, densities, viscosities, ix)
    nph = size(saturations, 1)
    if nph == 1
        rho = densities[ix]
        mu = viscosities[ix]
    else
        rho = zero(eltype(densities))
        mu = zero(eltype(viscosities))
        for ph in 1:nph
            s = saturations[ph, ix]
            rho += densities[ph, ix]*s
            mu += viscosities[ph, ix]*s
        end
    end
    return (rho, mu)
end

function Jutul.update_equation_in_entity!(eq_buf, self_cell, state, state0, eq::ConservationLaw{<:WellSegmentFlow}, model, dt, ldisc = local_discretization(eq, self_cell))
    (; cells, faces, signs) = ldisc

    mass = state.TotalMasses
    mass0 = state0.TotalMasses
    v = state.TotalMassFlux
    # For each component, compute fractional flow for mass flux + accumulation
    ncomp = number_of_components(model.system)
    m_t = component_sum(mass, self_cell)
    for i in 1:ncomp
        m_i = mass[i, self_cell]
        eq_i = (m_i - mass0[i, self_cell])/dt
        f_i = m_i/m_t
        for (cell, face, sgn) in zip(cells, faces, signs)
            v_f = sgn*v[face]
            f_o = mass[i, cell]/component_sum(mass, cell)
            eq_i += v_f*upw_flux(v_f, f_i, f_o)
        end
        eq_buf[i] = eq_i
    end
end

function component_sum(mass, i)
    s = zero(eltype(mass))
    for c = 1:size(mass, 1)
        s += mass[c, i]
    end
    return s
end
