
function single_unique_potential(model)
    # We should add capillary pressure here ventually
    return model.domain.discretizations.mass_flow.gravity
end

@inline function Jutul.compute_tpfa_flux!(q_i, left, right, face, face_sign, eq, state, model, dt, flow_disc)
    return component_mass_fluxes!(q_i, face, state, model, TPFA(left, right, face_sign), SPU(left, right))
end

@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:Union{ImmiscibleSystem, SinglePhaseSystem}, <:Any, <:Any}, kgrad, upw)
    for ph in eachindex(q)
        q_i = darcy_phase_mass_flux(face, ph, state, model, kgrad, upw)
        @inbounds q = @set q[ph] = q_i
    end
    return q
end

@inline function darcy_phase_mass_flux(face, phase, state, model, kgrad, upw)
    Q = darcy_phase_kgrad_potential(face, phase, state, model, kgrad)
    ρλ = c -> immiscible_phase_mass_mobility(state, phase, c)
    ρλ_f = upwind(upw, ρλ, Q)
    return ρλ_f*Q
end

@inline function darcy_phase_kgrad_potential(face, phase, state, model, tpfa::TPFA)
    trans = state.Transmissibilities
    grav = state.TwoPointGravityDifference
    ρ = state.PhaseMassDensities
    pc, ref_index = capillary_pressure(model, state)
    ∇p = pressure_gradient(state, tpfa)
    q = darcy_phase_flux_inner(∇p, pc, ref_index, trans, face, phase, grav, ρ, tpfa)
    return q
end

@inline function darcy_phase_flux_inner(∇p, pc, ref_index, trans, face, phase, grav, ρ, tpfa)
    l = tpfa.left
    r = tpfa.right
    s = tpfa.face_sign

    @inbounds T_f = trans[face]
    @inbounds gΔz = s*grav[face]

    Δpc = capillary_gradient(pc, l, r, phase, ref_index)

    @inbounds ρ_c = ρ[phase, l]
    @inbounds ρ_i = ρ[phase, r]
    ρ_avg = 0.5*(ρ_i + ρ_c)
    q = -T_f*(∇p + Δpc + gΔz*ρ_avg)
    return q
end

@inline function pressure_gradient(state, tpfa::TPFA)
    P = state.Pressure
    return @inbounds P[tpfa.left] - P[tpfa.right]
end

@inline function upwind(upw::SPU, F, q)
    if q >= 0
        up = upw.right
    else
        up = upw.left
    end
    return F(up)
end

@inline function immiscible_phase_mass_mobility(state, ph, c)
    @inbounds kr = state.RelativePermeabilities[ph, c]
    @inbounds ρ = state.PhaseMassDensities[ph, c]
    @inbounds μ = state.PhaseViscosities[ph, c]
    return ρ*kr/μ
end

capillary_gradient(::Nothing, c_l, c_r, ph, ph_ref) = 0.0
function capillary_gradient(pc, c_l, c_r, ph, ph_ref)
    if ph == ph_ref
        Δp_c = 0.0
    elseif ph < ph_ref
        Δp_c = pc[ph, c_l] - pc[ph, c_r]
    else
        Δp_c = pc[ph-1, c_l] - pc[ph-1, c_r]
    end
end



"""
TPFA KGrad(p) without gravity. (Outer version, with conn_data input)
"""
@inline function half_face_two_point_kgradp(conn_data::NamedTuple, p::AbstractArray)
    half_face_two_point_kgradp(conn_data.self, conn_data.other, conn_data.T, p)
end

"""
TPFA KGrad(p) without gravity. (Inner version, with explicit inputs)
"""
@inline function half_face_two_point_kgradp(c_self::I, c_other::I, T, p::AbstractArray{R}) where {R<:Real, I<:Integer}
    return -T*(p[c_self] - value(p[c_other]))
end

"""
TPFA-SPU Mobility * KGrad(p) without gravity. (Outer version, with conn_data input)
"""
@inline function half_face_two_point_flux_fused(conn_data, p, λ)
    return half_face_two_point_flux_fused(conn_data.self, conn_data.other, conn_data.T, p, λ)
end

"""
TPFA-SPU Mobility * KGrad(p) without gravity. (Inner version, with explicit inputs)
"""
@inline function half_face_two_point_flux_fused(c_self, c_other, T, p, λ)
    θ = half_face_two_point_kgradp(c_self, c_other, T, p)
    λᶠ = spu_upwind(c_self, c_other, θ, λ)
    return λᶠ*θ
end

"""
TPFA-SPU Mobility * (KGrad(p) + G). (Outer version, with conn_data input)
"""
@inline function half_face_two_point_flux_fused_gravity(conn_data, p, λ, density)
    return half_face_two_point_flux_fused_gravity(conn_data.self, conn_data.other, conn_data.T, p, λ, conn_data.gdz, density)
end

"""
TPFA-SPU Mobility * (KGrad(p) + G). (Inner version, with explicit inputs)
"""
@inline function half_face_two_point_flux_fused_gravity(c_self, c_other, T, p, λ, gΔz, density)
    θ = half_face_two_point_kgradp_gravity(c_self, c_other, T, p, gΔz, density)
    return spu_upwind_mult(c_self, c_other, θ, λ)
end


"""
Two point Darcy flux with gravity - outer version that takes in NamedTuple for static parameters
"""
@inline function half_face_two_point_kgradp_gravity(conn_data::NamedTuple, p, density)
    return half_face_two_point_kgradp_gravity(conn_data.self, conn_data.other, conn_data.T, p, conn_data.gdz, density)
end

"""
Two point Darcy flux with gravity - inner version that takes in cells and transmissibily explicitly
"""
@inline function half_face_two_point_kgradp_gravity(c_self::I, c_other::I, T, p::AbstractArray{R}, gΔz, ρ::AbstractArray{R}) where {R<:Real, I<:Integer}
    v = -T*Jutul.two_point_potential_drop_half_face(c_self, c_other, p, gΔz, ρ)
    return v
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
