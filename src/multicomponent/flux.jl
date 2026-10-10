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

@inline function inner_compositional!(q::SVector{M}, X, Y, q_l, q_v, upw, ::Val{N}) where {M, N}
    # Static component indices avoid SVector copies and large GPU local-memory buffers.
    return typeof(q)(ntuple(Val(M)) do i
        if i <= N
            X_f = upwind(upw, cell -> @inbounds(X[i, cell]), q_l)
            Y_f = upwind(upw, cell -> @inbounds(Y[i, cell]), q_v)
            q_l*X_f + q_v*Y_f
        else
            q[i]
        end
    end)
end

function add_diffusive_component_flux(q, ::Nothing, ::Nothing, face, state, model, grad, component_count)
    return q
end

function add_diffusive_component_flux(q, Dl, Dv, face, state, model, grad, component_count)
    sys = model.system
    masses = MultiComponentFlash.molar_masses(sys.equation_of_state)
    liquid = state.LiquidMassFractions
    vapor = state.VaporMassFractions
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

@inline function mass_to_mole_fractions(X, masses, cell, ::Val{N}) where N
    amounts = SVector{N}(ntuple(i -> @inbounds(X[i, cell]/masses[i]), Val(N)))
    return amounts/sum(amounts)
end

@inline diffusive_mole_fraction_accessor(X, masses, grad, count) =
    cell -> mass_to_mole_fractions(X, masses, cell, count)

@inline function diffusive_mole_fraction_accessor(X, masses, grad::TPFA, count)
    # Convert once per cell on a two-point face, rather than once per component.
    left, right = Jutul.cell_pair(grad)
    x_left = mass_to_mole_fractions(X, masses, left, count)
    x_right = mass_to_mole_fractions(X, masses, right, count)
    return cell -> if cell == left
        x_left
    else
        x_right
    end
end

@inline function add_phase_diffusive_component_flux(q::SVector{M}, D::AbstractMatrix, face, grad, S, ρ, mass_fractions::AbstractMatrix, phase, masses, ::Val{N}) where {M, N}
    fractions = diffusive_mole_fraction_accessor(mass_fractions, masses, grad, Val(N))
    molar_density = phase_diffusive_molar_density(ρ, phase, S, grad, fractions, masses)
    return typeof(q)(ntuple(Val(M)) do i
        @inbounds if i <= N
            fraction = cell -> fractions(cell)[i]
            molar_flux = -D[i, face]*molar_density*gradient(fraction, grad)
            q[i] + masses[i]*molar_flux
        else
            q[i]
        end
    end)
end

function phase_diffusive_molar_density(ρ, phase, S, grad, fractions, masses)
    density = cell -> @inbounds ρ[phase, cell]/dot(masses, fractions(cell))
    left, right = Jutul.cell_pair(grad)
    saturation = min(S[phase, left], S[phase, right])
    return saturation*face_average(density, grad)
end
