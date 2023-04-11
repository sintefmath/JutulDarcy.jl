@inline function Jutul.face_flux!(q_i, left, right, face, face_sign, eq::ConservationLaw{:TotalMasses, <:Any}, state, model::DarcyFlowModel, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    kgrad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    ft = Jutul.flux_type(eq)
    return component_mass_fluxes!(q_i, face, state, model, ft, kgrad, upw)
end

@inline function Jutul.face_flux!(q_i, face, eq::ConservationLaw{:TotalMasses, <:Any}, state, model::DarcyFlowModel, dt, flow_disc::PotentialFlow, ldisc)
    # Inner version, for generic flux
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    return component_mass_fluxes!(q_i, face, state, model, ft, kgrad, upw)
end

@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:Union{ImmiscibleSystem, SinglePhaseSystem}, <:Any, <:Any}, flux_type, kgrad, upw)
    disc = kgrad_common(face, state, model, kgrad)
    for ph in eachindex(q)
        q_i = darcy_phase_mass_flux(face, ph, state, model, kgrad, upw, disc)
        @inbounds q = setindex(q, q_i, ph)
    end
    return q
end

@inline function darcy_phase_mass_flux(face, phase, state, model, kgrad, upw, arg...)
    Q = darcy_phase_kgrad_potential(face, phase, state, model, kgrad, arg...)
    ρλ = state.PhaseMassMobilities
    ρλ_f = phase_upwind(upw, ρλ, phase, Q)
    return ρλ_f*Q
end

@inline function kgrad_common(face, state, model, tpfa::TPFA)
    ∇p = pressure_gradient(state, tpfa)
    trans = state.Transmissibilities
    grav = state.TwoPointGravityDifference
    @inbounds T_f = trans[face]
    @inbounds gΔz = tpfa.face_sign*grav[face]
    return (∇p, T_f, gΔz)
end

@inline function darcy_phase_kgrad_potential(face, phase, state, model, tpfa::TPFA, common = kgrad_common(face, state, model, tpfa))
    ρ = state.PhaseMassDensities
    pc, ref_index = capillary_pressure(model, state)
    ∇p, T_f, gΔz = common
    l = tpfa.left
    r = tpfa.right

    Δpc = capillary_gradient(pc, l, r, phase, ref_index)
    @inbounds ρ_c = ρ[phase, l]
    @inbounds ρ_i = ρ[phase, r]
    ρ_avg = 0.5*(ρ_i + ρ_c)
    q = -T_f*(∇p + Δpc + gΔz*ρ_avg)
    return q
end

@inline function gradient(X, tpfa::TPFA)
    return @inbounds X[tpfa.left] - X[tpfa.right]
end

pressure_gradient(state, disc) = gradient(state.Pressure, disc)

@inline function upwind(upw::SPU, F, q)
    flag = q >= 0
    if flag
        up = upw.right
    else
        up = upw.left
    end
    return F(up)
end

@inline function phase_upwind(upw, m::AbstractMatrix, phase::Integer, q)
    F(cell) = @inbounds m[phase, cell]
    return upwind(upw, F, q)
end

@inline capillary_gradient(::Nothing, c_l, c_r, ph, ph_ref) = 0.0
@inline function capillary_gradient(pc, c_l, c_r, ph, ph_ref)
    if ph == ph_ref
        Δp_c = zero(eltype(pc))
    else
        pos = ph - (ph > ph_ref)
        Δp_c = @inbounds pc[pos, c_l] - pc[pos, c_r]
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

@inline function saturation_averaged_density(ρ, ph, sat, c1, c2)
    @inbounds ρ_1 = ρ[ph, c1]
    @inbounds ρ_2 = ρ[ph, c2]
    @inbounds S_1 = sat[ph, c1]
    @inbounds S_2 = sat[ph, c2]

    avg = (ρ_1*S_1 + ρ_2*S_2)/max(S_1 + S_2, 1e-12)
    return avg
end

export component_mass_fluxes!, update_total_masses!
"""
    component_mass_fluxes!(q, face, state, model, kgrad, upw)

Implementation of component fluxes for a given system for a given face.
"""
function component_mass_fluxes!

end

"""
    update_total_masses!(totmass, tv, model, arg..., ix)

Update total masses for a given system. Number of input arguments varies based
on physical system under consideration.
"""
function update_total_masses!

end
