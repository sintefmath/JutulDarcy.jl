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

@inline function component_mass_fluxes!(q, face, state, model, flux_type, kgrad, upw)
    disc = flux_primitives(face, state, model, flux_type, kgrad, upw)
    for ph in eachindex(q)
        q_i = darcy_phase_mass_flux(face, ph, state, model, flux_type, kgrad, upw, disc)
        @inbounds q = setindex(q, q_i, ph)
    end
    return q
end

@inline function darcy_phase_mass_flux(face, phase, state, model, flux_type, kgrad, upw, arg...)
    Q = darcy_phase_kgrad_potential(face, phase, state, model, flux_type, kgrad, upw, arg...)
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

@inline function flux_primitives(face, state, model, flux_type::Jutul.DefaultFlux, tpfa::TPFA, upw)
    return kgrad_common(face, state, model, tpfa)
end

@inline function darcy_phase_kgrad_potential(face, phase, state, model, flux_type, tpfa::TPFA, upw, common = flux_primitives(face, state, model, flux_type, upw, tpfa))
    pc, ref_index = capillary_pressure(model, state)
    ∇p, T_f, gΔz = common
    l = tpfa.left
    r = tpfa.right

    Δpc = capillary_gradient(pc, l, r, phase, ref_index)
    ρ_avg = face_average_density(model, state, tpfa, phase)
    q = -T_f*(∇p + Δpc + gΔz*ρ_avg)
    return q
end

function face_average_density(model, state, tpfa, phase)
    ρ = state.PhaseMassDensities
    l = tpfa.left
    r = tpfa.right
    @inbounds ρ_c = ρ[phase, l]
    @inbounds ρ_i = ρ[phase, r]
    return 0.5*(ρ_i + ρ_c)
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
    return 0.5*(F(tpfa.right) + F(tpfa.left))
end

function phase_face_average(phase_property, tpfa, cell)
    F(cell) = @inbounds phase_property[phase, cell]
    return face_average(F, tpfa)
end

pressure_gradient(state, disc) = gradient(state.Pressure, disc)

@inline function upwind(upw::SPU, F, q)
    flag = q < 0
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
