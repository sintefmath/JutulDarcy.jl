"""
Half face Darcy flux (Blackoil version)
"""
function update_half_face_flux!(flux::AbstractArray, state, model::SimulationModel{D, S}, param, dt, flow_disc::TwoPointPotentialFlow{U, K, T}) where {D,S<:BlackOilSystem,U,K,T<:DarcyMassMobilityFlow}
    Rs = state.Rs
    kr = state.RelativePermeabilities
    μ = state.PhaseViscosities
    ρ = state.PhaseMassDensities
    P = state.Pressure
    b = state.ShrinkageFactors
    rhoS = tuple(param[:reference_densities]...)

    conn_data = flow_disc.conn_data
    pc, ref_index = capillary_pressure(model, state)

    sys = model.system
    blackoil_fluxes!(flux, conn_data, P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, model.context)
end

function blackoil_fluxes!(flux, conn_data, P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, context)
    tb = thread_batch(context)
    nf = size(flux, 2)
    # @batch minbatch = tb for i = 1:nf
    for i = 1:nf
        @inbounds @views blackoil_flux!(flux[:, i], conn_data[i], P, Rs, ρ, b, kr, μ, pc, ref_index, rhoS, (1,2,3))
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
