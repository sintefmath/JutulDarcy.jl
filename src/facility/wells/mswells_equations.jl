
"""
Hagedorn and Brown well bore friction model for a segment.
"""
struct SegmentWellBoreFrictionHB{R}
    L::R
    roughness::R
    D_outer::R
    D_inner::R
    assume_turbulent::Bool
    laminar_limit::R
    turbulent_limit::R
    function SegmentWellBoreFrictionHB(L, roughness, D_outer; D_inner = 0, assume_turbulent = false, laminar_limit = 2000.0, turbulent_limit = 4000.0)
        new{typeof(L)}(L, roughness, D_outer, D_inner, assume_turbulent, laminar_limit, turbulent_limit)
    end
end

function is_turbulent_flow(f::SegmentWellBoreFrictionHB, Re)
    return f.assume_turbulent || Re >= f.turbulent_limit
end

function is_laminar_flow(f::SegmentWellBoreFrictionHB, Re)
    return !f.assume_turbulent && Re <= f.laminar_limit
end

well_segment_is_closed(::SegmentWellBoreFrictionHB) = false

function segment_pressure_drop(f::SegmentWellBoreFrictionHB, v, ρ, μ)
    D⁰, Dⁱ = f.D_outer, f.D_inner
    R, L = f.roughness, f.L
    ΔD = D⁰-Dⁱ
    A = π*((D⁰/2)^2 - (Dⁱ/2)^2)
    # Scaling fix
    s = v > 0.0 ? 1.0 : -1.0
    e = eps(Float64)
    v = s*max(abs(v), e)

    Re = abs(D⁰*v/(A*μ))
    # Friction model - empirical relationship
    Re_l, Re_t = f.laminar_limit, f.turbulent_limit
    if is_laminar_flow(f, Re)
        f = 16.0/Re
    else
        # Either turbulent or intermediate flow regime. We need turbulent value either way.
        f_t = (-3.6*log10(6.9/Re +(R/(3.7*D⁰))^(10.0/9.0)))^(-2.0)
        if is_turbulent_flow(f, Re)
            # Turbulent flow
            f = f_t
        else
            # Intermediate regime - interpolation
            f_l = 16.0/Re_l
            Δf = f_t - f_l
            ΔRe = Re_t - Re_l
            f = f_l + (Δf / ΔRe)*(Re - Re_l)
        end
    end
    Δp = 2*f*L*v^2/((A^2)*D⁰*ρ)
    return Δp
end

struct ClosedSegment

end

well_segment_is_closed(::ClosedSegment) = true

"""
    PotentialDropBalanceWell(discretization)

Equation for the pressure drop equation in a multisegment well. This equation
lives on the segment between each node and balances the pressure difference
across the segment with the hydrostatic difference and well friction present in
the current flow regime.
"""
struct PotentialDropBalanceWell{T} <: JutulEquation
    flow_discretization::T
end

associated_entity(::PotentialDropBalanceWell) = Faces()

# local_discretization(e::PotentialDropBalanceWell, i) = e.flow_discretization(i, Faces())

Jutul.discretization(e::PotentialDropBalanceWell) = e.flow_discretization

import Jutul: two_point_potential_drop
function Jutul.update_equation_in_entity!(eq_buf, i, state, state0, eq::PotentialDropBalanceWell, model, dt, ldisc = local_discretization(eq, i))
    (; face, left, right, gdz) = ldisc
    μ = state.PhaseViscosities
    V = state.TotalMassFlux[face]
    densities = state.PhaseMassDensities
    s = state.Saturations
    p = state.Pressure

    w = physical_representation(model.domain)
    seg_model = w.segment_models[face]
    if well_segment_is_closed(seg_model)
        eq = V
    else
        rho_l, mu_l = saturation_mixed(s, densities, μ, left)
        rho_r, mu_r = saturation_mixed(s, densities, μ, right)
        rho = 0.5*(rho_l + rho_r)
        μ_mix = 0.5*(mu_l + mu_r)
        Δp = segment_pressure_drop(seg_model, V, rho, μ_mix)
        Δθ = two_point_potential_drop(p[left], p[right], gdz, rho_l, rho_r)
        eq = Δθ + Δp - 1e-12*V
    end
    eq_buf[] = eq
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
            if s > 0
                rho += densities[ph, ix]*s
                mu += viscosities[ph, ix]*s
            end
        end
    end
    return (rho, mu)
end

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::ConservationLaw{:TotalMasses, <:WellSegmentFlow}, model, dt, ldisc = local_discretization(eq, self_cell)) where T_e
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

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::ConservationLaw{:TotalThermalEnergy, <:WellSegmentFlow}, model, dt, ldisc = local_discretization(eq, self_cell)) where T_e
    (; cells, faces, signs) = ldisc
    nph = number_of_phases(model.system)
    λm = state.MaterialThermalConductivities
    mass = state.TotalMasses
    energy = state.TotalThermalEnergy
    energy0 = state0.TotalThermalEnergy
    density = state.PhaseMassDensities
    S = state.Saturations
    H_f = state.FluidEnthalpy
    v = state.TotalMassFlux
    T = state.Temperature

    eq = (energy[self_cell] - energy0[self_cell])/dt
    m_t_self = total_density(S, density, self_cell)

    for (cell, face, sgn) in zip(cells, faces, signs)
        v_f = sgn*v[face]
        for ph in 1:nph
            # Compute mass fraction of phase, then weight by enthalpy.
            m_t_other = total_density(S, density, cell)
            f_self = density[ph, self_cell]*S[ph, self_cell]/m_t_self
            f_other = density[ph, cell]*S[ph, cell]/m_t_other

            ME_other = H_f[ph, cell]*f_other
            ME_self = H_f[ph, self_cell]*f_self
            eq += v_f*upw_flux(v_f, ME_self, ME_other)
        end

        λm_f = λm[face]
        if λm_f >= 0.0
            # Account for heat conduction in well material
            eq += -λm_f*(T[cell] - T[self_cell])
        end
    end

    eq_buf[] = eq
end

function total_density(s, rho, cell)
    rho_tot = 0.0
    for ph in 1:size(rho, 1)
        rho_tot += rho[ph, cell]*s[ph, cell]
    end
    return rho_tot
end

function component_sum(mass, i)
    s = zero(eltype(mass))
    for c = axes(mass, 1)
        s += mass[c, i]
    end
    return s
end

function convergence_criterion(model, storage, eq::PotentialDropBalanceWell, eq_s, r; dt = 1.0, update_report = missing)
    e = (norm(r, Inf)/1e5, ) # Given as pressure - scale by 1 bar
    R = (AbsMax = (errors = e, names = "R"), )
    return R
end
