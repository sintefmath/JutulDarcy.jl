"""
Half face Darcy flux (Compositional version)
"""
function update_half_face_flux!(flux::AbstractArray, state, model::SimulationModel{D, S}, param, dt, flow_disc::TwoPointPotentialFlowHardCoded) where {D,S<:CompositionalSystem}
    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    kr = state.RelativePermeabilities
    μ = state.PhaseViscosities
    ρ = state.PhaseMassDensities
    P = state.Pressure
    Sat = state.Saturations

    conn_data = flow_disc.conn_data
    pc, ref_index = capillary_pressure(model, state)

    sys = model.system
    phase_ix = phase_indices(sys)
    aqua = Val(has_other_phase(sys))
    compositional_fluxes!(flux, conn_data, P, X, Y, ρ, kr, Sat, μ, pc, ref_index, aqua, phase_ix, model.context)
end

function compositional_fluxes!(flux, conn_data, P, X, Y, ρ, kr, Sat, μ, pc, ref_index, aqua, phase_ix, context)
    tb = minbatch(context)
    nf = size(flux, 2)
    @batch minbatch = tb for i = 1:nf
        @inbounds @views compositional_flux_gravity!(flux[:, i], conn_data[i], P, X, Y, ρ, kr, Sat, μ, pc, ref_index, aqua, phase_ix)
    end
end

function compositional_flux_gravity!(q, cd, P, X, Y, ρ, kr, Sat, μ, pc, ref_index, aqua, phase_ix)
    c, i, T = cd.self, cd.other, cd.T
    if haskey(cd, :gdz)
        gΔz = cd.gdz
    else
        gΔz = nothing
    end
    ∂ = (x) -> local_ad(x, c)
    return compute_compositional_flux_gravity!(q, c, i, ∂(P), ∂(X), ∂(Y), ∂(ρ), ∂(kr), ∂(Sat), ∂(μ), ∂(pc), T, gΔz, ref_index, aqua, phase_ix)
end

function compute_compositional_flux_gravity!(q, c, i, P, X, Y, ρ, kr, Sat, μ, pc, T, gΔz, ref_index, aqua::Val{false}, phase_ix)
    l, v = phase_ix

    G(phase) = compositional_buoyancy(ρ, phase, Sat, c, i, gΔz)
    Δpc(phase) = capillary_gradient(pc, c, i, phase, ref_index)
    mass_flux(Ψ, phase) = phase_mass_flux(Ψ, c, i, ρ, kr, μ, phase)

    @inbounds Δp = P[c] - P[i]

    ψ_l = -T*(Δp + Δpc(l) + G(l))
    ψ_v = -T*(Δp + Δpc(v) + G(v))

    F_l, c_l = mass_flux(ψ_l, l)
    F_v, c_v = mass_flux(ψ_v, v)

    for i in eachindex(q)
        @inbounds q[i] = F_l*X[i, c_l] + F_v*Y[i, c_v]
    end
end


function compute_compositional_flux_gravity!(q, c, i, P, X, Y, ρ, kr, Sat, μ, pc, T, gΔz, ref_index, aqua::Val{true}, phase_ix)
    a, l, v = phase_ix

    G(phase) = compositional_buoyancy(ρ, phase, Sat, c, i, gΔz)
    Δpc(phase) = capillary_gradient(pc, c, i, phase, ref_index)
    mass_flux(Ψ, phase) = phase_mass_flux(Ψ, c, i, ρ, kr, μ, phase)

    @inbounds Δp = P[c] - P[i]

    ψ_a = -T*(Δp + Δpc(a) + G(a))
    ψ_l = -T*(Δp + Δpc(l) + G(l))
    ψ_v = -T*(Δp + Δpc(v) + G(v))

    F_a, _ = mass_flux(ψ_a, a)
    F_l, c_l = mass_flux(ψ_l, l)
    F_v, c_v = mass_flux(ψ_v, v)

    N = length(q)
    for i in 1:(N-1)
        @inbounds q[i] = F_l*X[i, c_l] + F_v*Y[i, c_v]
    end
    @inbounds q[N] = F_a
end

compositional_buoyancy(ρ, phase, Sat, c, i, gΔz) = gΔz*saturation_averaged_density(ρ, phase, Sat, c, i)
compositional_buoyancy(ρ, phase, Sat, c, i, gΔz::Nothing) = zero(eltype(Sat))