@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:BlackOilSystem, <:Any, <:Any}, kgrad, upw)
    sys = model.system
    a, l, v = phase_indices(sys)
    rhoS = reference_densities(sys)
    # Aqueous phase is straightforward, just the mass flux for phase
    q_a = darcy_phase_mass_flux(face, a, state, model, kgrad, upw)
    # Get the potentials since the flux is more complicated for miscible phases
    ψ_l = darcy_phase_kgrad_potential(face, l, state, model, kgrad)
    ψ_v = darcy_phase_kgrad_potential(face, v, state, model, kgrad)

    kr = state.RelativePermeabilities
    μ = state.PhaseViscosities
    b = state.ShrinkageFactors

    # Oil
    f_bl = (cell, flag) -> b_mob(b, kr, μ, l, cell, flag)
    λb_l = upwind(upw, f_bl, ψ_l)
    q_l = rhoS[l]*λb_l*ψ_l
    # Gas
    f_bv = (cell, flag) -> b_mob(b, kr, μ, v, cell, flag)
    λb_v = upwind(upw, f_bv, ψ_v)

    f_rs = (cell, flag) -> @inbounds b[v, cell]
    # Note: Rs is upwinded by liquid potential
    rs = upwind(upw, f_rs, ψ_l)
    q_v = (λb_v*ψ_v + rs*λb_l*ψ_l)*rhoS[v]

    q = setindex(q, q_a, a)
    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end

function b_mob(b, kr, μ, ph, c, flag)
    @inbounds λ = kr[ph, c]/μ[ph, c]
    @inbounds b_f = b[ph, c]
    return λ*b_f
end

"""
Half face Darcy flux (Blackoil version)
"""
function update_half_face_flux!(flux::AbstractArray, state, model::SimulationModel{D, S}, param, dt, flow_disc::TwoPointPotentialFlowHardCoded) where {D,S<:BlackOilSystem}
    Rs = state.Rs
    kr = state.RelativePermeabilities
    μ = state.PhaseViscosities
    ρ = state.PhaseMassDensities
    P = state.Pressure
    b = state.ShrinkageFactors

    conn_data = flow_disc.conn_data
    pc, ref_index = capillary_pressure(model, state)

    sys = model.system
    ind = phase_indices(sys)
    rhoS = reference_densities(sys)
    blackoil_fluxes!(flux, conn_data, P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, ind, model.context)
end

function blackoil_fluxes!(flux, conn_data, P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, ind, context)
    tb = minbatch(context)
    nf = size(flux, 2)
    # @batch minbatch = tb for i = 1:nf
    for i = 1:nf
        @inbounds @views blackoil_flux!(flux[:, i], conn_data[i], P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, ind)
    end
end

function blackoil_flux!(q, cd, P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, phase_ix)
    c, i, T = cd.self, cd.other, cd.T
    if haskey(cd, :gdz)
        gΔz = cd.gdz
    else
        gΔz = nothing
    end
    ∂ = (x) -> local_ad(x, c)
    return blackoil_flux_internal!(q, c, i, ∂(P), ∂(Rs), ∂(ρ), ∂(b), ∂(kr), ∂(μ), ∂(pc), T, gΔz, ref_index, rhoS, phase_ix)
end

function blackoil_flux_internal!(q, c, i, P, Rs, ρ, b, kr, μ, pc, T, gΔz, ref_index, rhoS, phase_ix)
    a, l, v = phase_ix

    G(phase) = immiscible_buoyancy(ρ, phase, c, i, gΔz)
    Δpc(phase) = capillary_gradient(pc, c, i, phase, ref_index)
    # Note: b-factor weighted flux
    mass_flux(Ψ, phase) = phase_mass_flux(Ψ, c, i, b, kr, μ, phase)

    @inbounds Δp = P[c] - P[i]

    ψ_a = -T*(Δp + Δpc(a) + G(a))
    ψ_l = -T*(Δp + Δpc(l) + G(l))
    ψ_v = -T*(Δp + Δpc(v) + G(v))

    F_a, c_a = mass_flux(ψ_a, a)
    F_l, c_l = mass_flux(ψ_l, l)
    F_v, c_v = mass_flux(ψ_v, v)

    q[a] = F_a*rhoS[a]
    q[l] = F_l*rhoS[l]
    q[v] = (F_v + F_l*Rs[c_l])*rhoS[v]
end


immiscible_buoyancy(ρ, phase, c, i, gΔz) = gΔz*averaged_density(ρ, phase, c, i)
immiscible_buoyancy(ρ, phase, c, i, gΔz::Nothing) = zero(eltype(ρ))

@inline function averaged_density(ρ, ph, c1, c2)
    @inbounds ρ_1 = ρ[ph, c1]
    @inbounds ρ_2 = ρ[ph, c2]
    return (ρ_1 + ρ_2)/2
end
