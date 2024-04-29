"""
    PhaseMassFractions(:liquid)

Variable that defines the component mass fractions in a specific phase.
"""
struct PhaseMassFractions{T} <: CompositionalFractions
    phase::T
end

@jutul_secondary function update_phase_xy!(X, m::PhaseMassFractions, model::SimulationModel{D,S}, FlashResults, ix) where {D,S<:CompositionalSystem}
    eos = model.system.equation_of_state
    molar_mass = MultiComponentFlash.molar_masses(eos)
    n = MultiComponentFlash.number_of_components(eos)
    phase = m.phase
    @inbounds for i in ix
        f = FlashResults[i]
        if phase_is_present(phase, f.state)
            # X_i = view(X, :, i)
            r = phase_data(f, phase)
            x_i = r.mole_fractions
            update_mass_fractions!(X, x_i, i, molar_mass, n)
        end
    end
end

@inline function update_mass_fractions!(X, x, cell, molar_masses, n)
    t = zero(eltype(X))
    @inbounds for i in 1:n
        t += molar_masses[i] * x[i]
    end
    @inbounds for i in 1:n
        X[i, cell] = molar_masses[i] * x[i] / t
    end
end

# Total masses
@jutul_secondary function update_total_masses!(
        totmass,
        tmvar::TotalMasses,
        model::LVCompositionalModel2Phase,
        FlashResults,
        PhaseMassDensities,
        Saturations,
        VaporMassFractions,
        LiquidMassFractions,
        FluidVolume,
        ix
    )
    compositional_mass_update_loop!(
        totmass,
        model,
        FlashResults,
        PhaseMassDensities,
        Saturations,
        LiquidMassFractions,
        VaporMassFractions,
        nothing,
        FluidVolume,
        ix
    )
end

@jutul_secondary function update_total_masses!(
        totmass,
        tmvar::TotalMasses,
        model::LVCompositionalModel3Phase,
        FlashResults,
        PhaseMassDensities,
        Saturations,
        VaporMassFractions,
        LiquidMassFractions,
        ImmiscibleSaturation,
        FluidVolume,
        ix
    )
    compositional_mass_update_loop!(
        totmass,
        model,
        FlashResults,
        PhaseMassDensities,
        Saturations,
        LiquidMassFractions,
        VaporMassFractions,
        ImmiscibleSaturation,
        FluidVolume,
        ix
    )
end


function compositional_mass_update_loop!(totmass, model, F, ρ, Sat, X, Y, sw, pv, ix)
    sys = model.system
    phase_ix = phase_indices(sys)
    N = size(totmass, 1)
    for cell in ix
        @inbounds two_phase_compositional_mass!(totmass, F[cell].state, sw, pv, ρ, X, Y, Sat, cell, N, phase_ix)
    end
    return totmass
end

function degrees_of_freedom_per_entity(model::SimulationModel{G,S}, v::TotalMasses) where {G<:Any,S<:MultiComponentSystem}
    number_of_components(model.system)
end

"""
Update total masses for two-phase compositional
"""
function two_phase_compositional_mass!(M, state, sw::Nothing, Φ, ρ, X, Y, S, cell, N, phase_ix)
    update_mass_two_phase_compositional!(M, state, sw, Φ, ρ, X, Y, S, cell, phase_ix, N)
end

"""
Update total masses for two-phase compositional where another immiscible phase is present
"""
function two_phase_compositional_mass!(M, state, sw, Φ, ρ, X, Y, S, cell, N, phase_ix)
    update_mass_two_phase_compositional!(M, state, sw, Φ, ρ, X, Y, S, cell, phase_ix[2:end], N - 1)
    a, = phase_ix
    @inbounds M[end, cell] = ρ[a, cell]*S[a, cell]*Φ[cell]
end

function update_mass_two_phase_compositional!(M, state, sw, Φ, ρ, X, Y, S, cell, phase_ix, N)
    l, v = phase_ix
    has_liquid = liquid_phase_present(state)
    has_vapor = vapor_phase_present(state)
    if has_liquid && has_vapor
        two_phase_mass!(M, ρ, S, X, Y, Φ, cell, N, l, v)
    elseif has_liquid
        single_phase_mass!(M, ρ, S, X, Φ, cell, N, l)
    else
        single_phase_mass!(M, ρ, S, Y, Φ, cell, N, v)
    end
end

function single_phase_mass!(M, ρ, S, mass_fractions, Φ, cell, N, phase)
    S_eos = S[phase, cell]
    if S_eos < MINIMUM_COMPOSITIONAL_SATURATION
        S_eos = replace_value(S_eos, MINIMUM_COMPOSITIONAL_SATURATION)
    end
    @inbounds M_l = ρ[phase, cell] * S_eos
    for c in 1:N
        @inbounds M[c, cell] = M_l*mass_fractions[c, cell]*Φ[cell]
    end
end

function two_phase_mass!(M, ρ, S, X, Y, Φ, cell, N, l, v)
    @inbounds M_l = ρ[l, cell] * S[l, cell]
    @inbounds M_v = ρ[v, cell] * S[v, cell]
    for c in 1:N
        @inbounds M[c, cell] = (M_l*X[c, cell] + M_v*Y[c, cell])*Φ[cell]
    end
end
