@inline function Jutul.face_flux!(q_i, left, right, face, face_sign, eq::ConservationLaw{:TotalMasses, <:Any}, state, model::DarcyFlowModel, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    kgrad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    ft = Jutul.flux_type(eq)
    return @inbounds component_mass_fluxes!(q_i, face, state, model, ft, kgrad, upw)
end

@inline function Jutul.face_flux!(q_i, face, eq::ConservationLaw{:TotalMasses, <:Any}, state, model::DarcyFlowModel, dt, flow_disc::PotentialFlow, ldisc)
    # Inner version, for generic flux
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    return @inbounds component_mass_fluxes!(q_i, face, state, model, ft, kgrad, upw)
end

@inline function component_mass_fluxes!(q, face, state, model, flux_type, kgrad, upw)
    T = eltype(q)
    phase_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, kgrad, upw)
    for ph in eachindex(phase_fluxes)
        q_ph = phase_fluxes[ph]
        @inbounds q = setindex(q, q_ph, ph)
    end
    return q
end

@inline function darcy_phase_mass_fluxes(face, state, model, flux_type, kgrad, upw, phases = eachphase(model.system))
    dpot = darcy_permeability_potential_differences(face, state, model, flux_type, kgrad, upw, phases)
    F(phase) = darcy_phase_mass_flux(face, phase, state, model, flux_type, kgrad, upw, dpot[phase])
    return map(
        F,
        phases
    )
end

@inline function darcy_phase_mass_flux(face, phase, state, model, flux_type, kgrad, upw, Q = missing)
    if ismissing(Q)
        Q = darcy_permeability_potential_differences(face, state, model, flux_type, kgrad, upw, phase)[1]
    end
    if haskey(state, :PhaseMassMobilities)
        ρλ = state.PhaseMassMobilities
        ρλ_f = phase_upwind(upw, ρλ, phase, Q)
    else
        rho = state.PhaseMassDensities
        λ = state.PhaseMobilities
        F = cell -> @inbounds λ[phase, cell]*rho[phase, cell]
        ρλ_f = upwind(upw, F, Q)
    end
    return ρλ_f*Q
end

@inline function darcy_phase_volume_fluxes(face, state, model, flux_type, kgrad, upw, phases = eachphase(model.system))
    dpot = darcy_permeability_potential_differences(face, state, model, flux_type, kgrad, upw, phases)
    F(phase) = darcy_phase_volume_flux(face, phase, state, model, flux_type, kgrad, upw, dpot[phase])
    return map(
        F,
        phases
    )
end

@inline function darcy_phase_volume_flux(face, phase, state, model, flux_type, kgrad, upw, Q = missing)
    if ismissing(Q)
        Q = darcy_phase_kgrad_potential(face, phase, state, model, flux_type, kgrad, upw, phase)
    end
    λ = state.PhaseMobilities
    F = cell -> @inbounds λ[phase, cell]
    ρλ_f = upwind(upw, F, Q)
    return ρλ_f*Q
end

@inline function darcy_permeability_potential_differences(
        face,
        state,
        model,
        flux_type,
        kgrad,
        upw,
        phases = eachphase(model.system)
    )
    T_f = effective_transmissibility(state, face, kgrad)
    gΔz = effective_gravity_difference(state, face, kgrad)
    pc, ref_index = capillary_pressure(model, state)

    ∇p = pressure_gradient(state, kgrad)

    @inline function phase_pot(phase)
        Δpc = capillary_gradient(pc, kgrad, phase, ref_index)
        ρ_avg = face_average_density(model, state, kgrad, phase)
        return -T_f*(∇p + Δpc + gΔz*ρ_avg)
    end
    return map(phase_pot, phases)
end

@inline function effective_transmissibility(state, face, kgrad)
    @inbounds T_f = state.Transmissibilities[face]
    if haskey(state, :PermeabilityMultiplier)
        K_mul = state.PermeabilityMultiplier
        get_kval(c) = @inbounds K_mul[c]
        m = face_average(get_kval, kgrad)
        T_f *= m
    end
    return T_f
end

function effective_gravity_difference(state, face, kgrad)
    grav = state.TwoPointGravityDifference
    face_sign(::Any) = 1
    face_sign(x::TPFA) = x.face_sign
    @inbounds gΔz = face_sign(kgrad)*grav[face]
    return gΔz
end

@inline function darcy_permeability_potential_difference(
        face,
        state,
        model,
        flux_type,
        kgrad,
        upw,
        phase::Int
    )
    dpots = darcy_permeability_potential_differences(
        face,
        state,
        model,
        flux_type,
        kgrad,
        upw,
        (phase, )
    )
    return dpots[1]
end

@inline function face_average_density(model, state, tpfa, phase, ρ = state.PhaseMassDensities)
    return phase_face_average(ρ, tpfa, phase)
end

@inline function face_average_density(model::CompositionalModel, state, tpfa, phase, ρ = state.PhaseMassDensities)
    sys = flow_system(model.system)
    l = tpfa.left
    r = tpfa.right
    @inbounds ρ_l = ρ[phase, l]
    @inbounds ρ_r = ρ[phase, r]

    s = state.Saturations
    ϵ = MINIMUM_COMPOSITIONAL_SATURATION
    @inbounds s_l = s[phase, l]
    @inbounds s_r = s[phase, r]

    s_l_tiny = s_l <= ϵ
    s_r_tiny = s_r <= ϵ
    if s_l_tiny && s_r_tiny
        ρ_avg = zero(s_l)
    elseif s_l_tiny
        ρ_avg = ρ_r
    elseif s_r_tiny
        ρ_avg = ρ_l
    else
        # alt def: (s_l*ρ_r + s_r*ρ_l)/(s_l + s_r)
        ρ_avg = 0.5*(ρ_r + ρ_l)
    end
    return ρ_avg
end

@inline function gradient(X::AbstractVector, tpfa::TPFA)
    return @inbounds X[tpfa.right] - X[tpfa.left]
end

@inline function gradient(X::AbstractMatrix, i, tpfa::TPFA)
    return @inbounds X[i, tpfa.right] - X[i, tpfa.left]
end

@inline function gradient(F, tpfa::TPFA)
    # Function handle version
    return F(tpfa.right) - F(tpfa.left)
end

function face_average(F, tpfa)
    l, r = Jutul.cell_pair(tpfa)
    return 0.5*(F(r) + F(l))
end

function phase_face_average(phase_property, tpfa, phase)
    F(cell) = @inbounds phase_property[phase, cell]
    return face_average(F, tpfa)
end

pressure_gradient(state, disc) = gradient(state.Pressure, disc)

@inline function upwind(upw::SPU, F, q)
    flag = q < zero(q)
    if flag
        up = upw.right
    else
        up = upw.left
    end
    return F(up)
end

@inline function upwind(upw::SPU, X::AbstractArray, q)
    flag = q < zero(q)
    if flag
        up = upw.right
    else
        up = upw.left
    end
    return @inbounds X[up]
end

@inline function phase_upwind(upw, m::AbstractMatrix, phase::Integer, q)
    F(cell) = @inbounds m[phase, cell]
    return upwind(upw, F, q)
end

@inline function upwind(upw::Jutul.WENO.WENOFaceDiscretization, F, q)
    return Jutul.WENO.weno_upwind(upw, F, q)
end

@inline function capillary_gradient(pc, tpfa::TPFA, ph, ph_ref)
    return capillary_gradient(pc, tpfa.left, tpfa.right, ph, ph_ref)
end

@inline function capillary_gradient(::Nothing, c_l, c_r, ph, ph_ref)
    return 0.0
end

@inline function capillary_gradient(pc, c_l, c_r, ph, ph_ref)
    if ph == ph_ref
        Δp_c = zero(eltype(pc))
    else
        pos = ph - (ph > ph_ref)
        Δp_c = @inbounds pc[pos, c_r] - pc[pos, c_l]
    end
    return Δp_c
end



@inline function phase_mass_flux(Ψ, c, i, ρ, kr, μ, ph)
    upc = upwind_cell(Ψ, c, i)
    @inbounds F = ρ[ph, upc]*(kr[ph, upc]/μ[ph, upc])*Ψ
    return (F, upc)
end

@inline function upwind_cell(pot, l, r)
    if pot < 0
        c = l
    else
        c = r
    end
end

"""
    component_mass_fluxes!(q, face, state, model, flux_type, kgrad, upw)

Implementation of component fluxes for a given system for a given face. Should
return a `StaticVector` with one entry per component.

$(SIGNATURES)

"""
function component_mass_fluxes!

end

"""
    update_total_masses!(totmass, tv, model, arg..., ix)

Update total masses for a given system. Number of input arguments varies based
on physical system under consideration.

$(SIGNATURES)

"""
function update_total_masses!

end
