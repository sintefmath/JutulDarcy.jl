@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:CompositionalSystem, <:Any, <:Any}, kgrad, upw)
    # for ph in eachindex(q)
    #     q_i = darcy_phase_mass_flux(face, ph, state, model, kgrad, upw)
    #     @inbounds q = setindex(q, q_i, ph)
    # end
    sys = model.system
    aqua = Val(has_other_phase(sys))
    ph_ix = phase_indices(sys)

    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    kr = state.RelativePermeabilities
    μ = state.PhaseViscosities
    ρ = state.PhaseMassDensities
    pressure = state.Pressure
    s = state.Saturations

    q = compositional_fluxes!(q, face, state,  X, Y, kr, μ, ρ, pressure, s, model, kgrad, upw, aqua, ph_ix)
    return q
end

@inline function compositional_fluxes!(q, face, state, X, Y, kr, μ, ρ, pressure, s, model, kgrad, upw, aqua::Val{false}, phase_ix)
    nc = size(X, 1)
    l, v = phase_ix
    q_l = darcy_phase_mass_flux(face, l, state, model, kgrad, upw)
    q_v = darcy_phase_mass_flux(face, v, state, model, kgrad, upw)

    q = inner_compositional!(q, X, Y, q_l, q_v, upw, nc)
    return q
end

@inline function compositional_fluxes!(q, face, state, X, Y, kr, μ, ρ, pressure, s, model, kgrad, upw, aqua::Val{true}, phase_ix)
    nc = size(X, 1)
    a, l, v = phase_ix
    q_a = darcy_phase_mass_flux(face, a, state, model, kgrad, upw)
    q_l = darcy_phase_mass_flux(face, l, state, model, kgrad, upw)
    q_v = darcy_phase_mass_flux(face, v, state, model, kgrad, upw)

    q = inner_compositional!(q, X, Y, q_l, q_v, upw, nc)
    q = setindex(q, q_a, nc+1)
    return q
end

@inline function inner_compositional!(q, X, Y, q_l, q_v, upw, nc)
    for i in 1:nc
        X_f = upwind(upw, cell -> @inbounds(X[i, cell]), q_l)
        Y_f = upwind(upw, cell -> @inbounds(Y[i, cell]), q_v)

        q_i = q_l*X_f + q_v*Y_f
        q = setindex(q, q_i, i)
    end
    return q
end

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