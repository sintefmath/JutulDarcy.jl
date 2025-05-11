@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:CompositionalSystem, <:Any, <:Any}, flux_type, kgrad, upw)
    sys = model.system
    aqua = Val(has_other_phase(sys))
    ph_ix = phase_indices(sys)

    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    S = state.Saturations
    ρ = state.PhaseMassDensities
    if haskey(state, :Diffusivities)
        D = state.Diffusivities
    else
        D = nothing
    end
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, kgrad, upw)
    q = compositional_fluxes!(q, face, state, S, ρ, X, Y, D, model, flux_type, kgrad, mass_fluxes, upw, aqua, ph_ix)
    return q
end

@inline function compositional_fluxes!(q, face, state, S, ρ, X, Y, D, model, flux_type, kgrad, mass_fluxes, upw, aqua::Val{false}, phase_ix)
    nc = size(X, 1)
    l, v = phase_ix
    q_l = mass_fluxes[l]
    q_v = mass_fluxes[v]

    q = inner_compositional!(q, S, ρ, X, Y, D, q_l, q_v, face, kgrad, upw, nc, phase_ix)
    return q
end

@inline function compositional_fluxes!(q, face, state, S, ρ, X, Y, D, model, flux_type, kgrad, mass_fluxes, upw, aqua::Val{true}, phase_ix)
    nc = size(X, 1)
    a, l, v = phase_ix
    q_a = mass_fluxes[a]
    q_l = mass_fluxes[l]
    q_v = mass_fluxes[v]

    q = inner_compositional!(q, S, ρ, X, Y, D, q_l, q_v, face, kgrad, upw, nc, (l, v))
    q = setindex(q, q_a, nc+1)
    return q
end

@inline function inner_compositional!(q, S, ρ, X, Y, D, q_l, q_v, face, grad, upw, nc, lv)
    for i in 1:nc
        X_f = upwind(upw, cell -> @inbounds(X[i, cell]), q_l)
        Y_f = upwind(upw, cell -> @inbounds(Y[i, cell]), q_v)
        q_i = q_l*X_f + q_v*Y_f
        q = setindex(q, q_i, i)
    end
    q = add_diffusive_component_flux(q, S, ρ, X, Y, D, face, grad, lv)
    return q
end

function add_diffusive_component_flux(q, S, ρ, X, Y, D::Nothing, face, grad, lv)
    return q
end

function add_diffusive_component_flux(q, S, ρ, X, Y, D, face, grad, lv)
    l, v = lv
    diff_mass_l = phase_diffused_mass(D, ρ, l, S, face, grad)
    diff_mass_v = phase_diffused_mass(D, ρ, v, S, face, grad)

    @inbounds for i in eachindex(q)
        q_i = q[i]
        dX = gradient(X, i, grad)
        dY = gradient(Y, i, grad)

        q_i += diff_mass_l*dX + diff_mass_v*dY
        q = setindex(q, q_i, i)
    end
    return q
end

function phase_diffused_mass(D, ρ, α, S, face, grad)
    @inbounds D_α = D[α, face]
    den_α = cell -> @inbounds ρ[α, cell]
    # Take minimum - diffusion should not cross phase boundaries.
    left, right = Jutul.cell_pair(grad)
    S = min(S[α, left], S[α, right])
    return -D_α*S*face_average(den_α, grad)
end
