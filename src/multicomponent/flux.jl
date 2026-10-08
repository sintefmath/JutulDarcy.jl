@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:CompositionalSystem, <:Any, <:Any}, flux_type, kgrad, upw)
    sys = model.system
    aqua = Val(has_other_phase(sys))
    component_count = Val(MultiComponentFlash.number_of_components(
        sys.equation_of_state))
    ph_ix = phase_indices(sys)

    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    S = state.Saturations
    ρ = state.PhaseMassDensities
    D = (phase_diffusivities(state, LiquidPhase()), phase_diffusivities(state, VaporPhase()))
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, kgrad, upw)
    q = compositional_fluxes!(q, face, state, S, ρ, X, Y, D, model,
        flux_type, kgrad, mass_fluxes, upw, aqua, ph_ix, component_count)
    return q
end

@inline function compositional_fluxes!(q, face, state, S, ρ, X, Y, D,
        model, flux_type, kgrad, mass_fluxes, upw, aqua::Val{false},
        phase_ix, component_count::Val{N}) where N
    l, v = phase_ix
    q_l = mass_fluxes[l]
    q_v = mass_fluxes[v]

    q = inner_compositional!(q, S, ρ, X, Y, D, q_l, q_v, face,
        kgrad, upw, component_count, phase_ix)
    return q
end

@inline function compositional_fluxes!(q, face, state, S, ρ, X, Y, D,
        model, flux_type, kgrad, mass_fluxes, upw, aqua::Val{true},
        phase_ix, component_count::Val{N}) where N
    a, l, v = phase_ix
    q_a = mass_fluxes[a]
    q_l = mass_fluxes[l]
    q_v = mass_fluxes[v]

    q = inner_compositional!(q, S, ρ, X, Y, D, q_l, q_v, face,
        kgrad, upw, component_count, (l, v))
    q = setindex(q, q_a, N + 1)
    return q
end

@inline function inner_compositional!(q, S, ρ, X, Y, D, q_l, q_v,
        face, grad, upw, ::Val{N}, lv) where N
    for i in 1:N
        X_f = upwind(upw, cell -> @inbounds(X[i, cell]), q_l)
        Y_f = upwind(upw, cell -> @inbounds(Y[i, cell]), q_v)
        q_i = q_l*X_f + q_v*Y_f
        q = setindex(q, q_i, i)
    end
    q = add_diffusive_component_flux(q, S, ρ, X, Y, D, face, grad,
        lv, Val(N))
    return q
end

function add_diffusive_component_flux(q, S, ρ, X, Y, D, face, grad, lv,
        ::Val{N}) where N
    l, v = lv
    Dl, Dv = D
    @inbounds for i in 1:N
        q_i = q[i]
        if !isnothing(Dl)
            diff_mass_l = phase_diffused_mass(Dl, ρ, l, i, S, face, grad)
            q_i += diff_mass_l*gradient(X, i, grad)
        end
        if !isnothing(Dv)
            diff_mass_v = phase_diffused_mass(Dv, ρ, v, i, S, face, grad)
            q_i += diff_mass_v*gradient(Y, i, grad)
        end
        q = setindex(q, q_i, i)
    end
    return q
end

function phase_diffused_mass(D, ρ, phase, component, S, face, grad)
    @inbounds coefficient = D[component, face]
    density = cell -> @inbounds ρ[phase, cell]
    # Take minimum - diffusion should not cross phase boundaries.
    left, right = Jutul.cell_pair(grad)
    saturation = min(S[phase, left], S[phase, right])
    return -coefficient*saturation*face_average(density, grad)
end
