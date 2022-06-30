# local_discretization(e::PotentialDropBalanceWell, i) = e.flow_discretization(i, Faces())

Jutul.discretization(e::PotentialDropBalanceWell) = e.flow_discretization

function Jutul.update_equation_in_entity!(eq_buf, i, state, state0, eq::PotentialDropBalanceWell, model, dt, ldisc = local_discretization(eq, i))
    (; face, left, right, gdz) = ldisc
    μ = state.PhaseViscosities
    V = state.TotalMassFlux[face]
    densities = state.PhaseMassDensities
    s = state.Saturations

    rho_l, mu_l = saturation_mixed(s, densities, μ, left)
    rho_r, mu_r = saturation_mixed(s, densities, μ, right)

    Δθ = Jutul.two_point_potential_drop(left, right, gdz, rho_l, rho_r)
    if Δθ > 0
        μ_mix = mu_l
    else
        μ_mix = mu_r
    end
    rho = 0.5*(rho_l + rho_r)

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

function Jutul.update_equation_in_entity!(eq_buf, i, state, state0, eq::ConservationLaw, model, dt, ldisc = local_discretization(eq, i))
    error()
    (; face, left, right, gdz) = ldisc
    μ = state.PhaseViscosities
    V = state.TotalMassFlux[face]
    densities = state.PhaseMassDensities
    s = state.Saturations

    rho_l, mu_l = saturation_mixed(s, densities, μ, left)
    rho_r, mu_r = saturation_mixed(s, densities, μ, right)

    Δθ = Jutul.two_point_potential_drop(left, right, gdz, rho_l, rho_r)
    if Δθ > 0
        μ_mix = mu_l
    else
        μ_mix = mu_r
    end
    rho = 0.5*(rho_l + rho_r)

    seg_model = model.domain.grid.segment_models[face]
    Δp = segment_pressure_drop(seg_model, V, rho, μ_mix)

    eq_buf[] = Δθ + Δp
end
