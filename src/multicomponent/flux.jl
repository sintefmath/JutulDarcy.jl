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
    if haskey(state, :Diffusivities)
        D = state.Diffusivities
    else
        D = nothing
    end
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, kgrad, upw)
    q = compositional_fluxes!(q, face, state, S, ρ, X, Y, D, model,
        flux_type, kgrad, mass_fluxes, upw, aqua, ph_ix, component_count)
    if haskey(state, :MolarDiffusivities)
        q = add_molar_diffusive_flux(q, face, state, model, kgrad, component_count)
    end
    return q
end

function add_molar_diffusive_flux(q, face, state, model, grad, ::Val{N}) where N
    eos = model.system.equation_of_state
    mw = MultiComponentFlash.molar_masses(eos)
    D = state.MolarDiffusivities
    left, right = Jutul.cell_pair(grad)
    for (phase, phase_name) in ((liquid_phase_index(model.system), Val(:liquid)),
                              (vapor_phase_index(model.system), Val(:vapor)))
        function concentration(cell)
            fractions = phase_data(state.FlashResults[cell], phase_name).mole_fractions
            return state.PhaseMassDensities[phase, cell]/sum(mw .* fractions)
        end
        c = face_average(concentration, grad)
        s = min(state.Saturations[phase, left], state.Saturations[phase, right])
        for i in 1:N
            mole_fraction(cell) = phase_data(state.FlashResults[cell], phase_name).mole_fractions[i]
            dq = -D[(phase - 1)*N + i, face]*s*c*mw[i]*gradient(mole_fraction, grad)
            q = setindex(q, q[i] + dq, i)
        end
    end
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

function add_diffusive_component_flux(q, S, ρ, X, Y, D::Nothing, face,
        grad, lv, component_count::Val)
    return q
end

function add_diffusive_component_flux(q, S, ρ, X, Y, D, face, grad, lv,
        ::Val{N}) where N
    l, v = lv
    diff_mass_l = phase_diffused_mass(D, ρ, l, S, face, grad)
    diff_mass_v = phase_diffused_mass(D, ρ, v, S, face, grad)

    @inbounds for i in 1:N
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
