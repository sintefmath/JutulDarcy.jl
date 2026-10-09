@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:CompositionalSystem, <:Any, <:Any}, flux_type, kgrad, upw)
    sys = model.system
    aqua = Val(has_other_phase(sys))
    component_count = Val(MultiComponentFlash.number_of_components(sys.equation_of_state))
    ph_ix = phase_indices(sys)

    X = state.LiquidMassFractions
    Y = state.VaporMassFractions
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, kgrad, upw)
    q = compositional_fluxes!(q, X, Y, mass_fluxes, upw, aqua, ph_ix, component_count)

    Dl = phase_diffusivities(state, LiquidPhase())
    Dv = phase_diffusivities(state, VaporPhase())
    return add_diffusive_component_flux(q, Dl, Dv, face, state, model, kgrad, component_count)
end

@inline function compositional_fluxes!(q, X, Y, mass_fluxes, upw, ::Val{false},
        phase_ix, component_count)
    l, v = phase_ix
    return inner_compositional!(q, X, Y, mass_fluxes[l], mass_fluxes[v], upw, component_count)
end

@inline function compositional_fluxes!(q, X, Y, mass_fluxes, upw, ::Val{true},
        phase_ix, component_count::Val{N}) where N
    a, l, v = phase_ix
    q = inner_compositional!(q, X, Y, mass_fluxes[l], mass_fluxes[v], upw, component_count)
    return setindex(q, mass_fluxes[a], N + 1)
end

@inline function inner_compositional!(q, X, Y, q_l, q_v, upw, ::Val{N}) where N
    for i in 1:N
        X_f = upwind(upw, cell -> @inbounds(X[i, cell]), q_l)
        Y_f = upwind(upw, cell -> @inbounds(Y[i, cell]), q_v)
        q_i = q_l*X_f + q_v*Y_f
        q = setindex(q, q_i, i)
    end
    return q
end

function add_diffusive_component_flux(q, ::Nothing, ::Nothing, face, state, model, grad, component_count)
    return q
end

function add_diffusive_component_flux(q, Dl, Dv, face, state, model, grad, component_count)
    sys = model.system
    masses = MultiComponentFlash.molar_masses(sys.equation_of_state)
    liquid = cell -> phase_data(state.FlashResults[cell], Val(:liquid)).mole_fractions
    vapor = cell -> phase_data(state.FlashResults[cell], Val(:vapor)).mole_fractions
    S = state.Saturations
    ρ = state.PhaseMassDensities
    q = add_phase_diffusive_component_flux(q, Dl, face, grad, S, ρ,
        liquid, liquid_phase_index(sys), masses, component_count)
    q = add_phase_diffusive_component_flux(q, Dv, face, grad, S, ρ,
        vapor, vapor_phase_index(sys), masses, component_count)
    return q
end

function add_phase_diffusive_component_flux(q, ::Nothing, face, grad, S, ρ, fractions, phase, masses, component_count::Val)
    return q
end

function add_phase_diffusive_component_flux(q, D, face, grad, S, ρ, fractions, phase, masses, ::Val{N}) where N
    molar_density = phase_diffusive_molar_density(ρ, phase, S, grad, fractions, masses)
    @inbounds for i in 1:N
        fraction = cell -> fractions(cell)[i]
        molar_flux = -D[i, face]*molar_density*gradient(fraction, grad)
        q = setindex(q, q[i] + masses[i]*molar_flux, i)
    end
    return q
end

function phase_diffusive_molar_density(ρ, phase, S, grad, fractions, masses)
    density = cell -> @inbounds ρ[phase, cell]/dot(masses, fractions(cell))
    left, right = Jutul.cell_pair(grad)
    saturation = min(S[phase, left], S[phase, right])
    return saturation*face_average(density, grad)
end
