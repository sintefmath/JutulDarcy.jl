struct PhaseMassFractions{T} <: CompositionalFractions
    phase::T
end

@jutul_secondary function update_phase_xy!(X, m::PhaseMassFractions, model::SimulationModel{D,S}, FlashResults, ix) where {D,S<:CompositionalSystem}
    molar_mass = map((x) -> x.mw, model.system.equation_of_state.mixture.properties)
    phase = m.phase
    @inbounds for i in ix
        f = FlashResults[i]
        if phase_is_present(phase, f.state)
            # X_i = view(X, :, i)
            r = phase_data(f, phase)
            x_i = r.mole_fractions
            update_mass_fractions!(X, x_i, i, molar_mass)
        end
    end
end

@inline function update_mass_fractions!(X, x, cell, molar_masses)
    t = zero(eltype(X))
    @inbounds for i in eachindex(x)
        t += molar_masses[i] * x[i]
    end
    @inbounds for i in eachindex(x)
        X[i, cell] = molar_masses[i] * x[i] / t
    end
end

# Total masses
@jutul_secondary function update_total_masses!(totmass, tv::TotalMasses, model::SimulationModel{G,S},
                                                                                                    FlashResults,
                                                                                                    PhaseMassDensities,
                                                                                                    Saturations,
                                                                                                    VaporMassFractions,
                                                                                                    LiquidMassFractions,
                                                                                                    FluidVolume, ix) where {G,S<:CompositionalSystem}
    pv = FluidVolume
    ρ = PhaseMassDensities
    X = LiquidMassFractions
    Y = VaporMassFractions
    Sat = Saturations
    F = FlashResults
    sys = model.system
    phase_ix = phase_indices(sys)
    has_other = Val(has_other_phase(sys))
    N = size(totmass, 1)
    for cell in ix
        @inbounds two_phase_compositional_mass!(totmass, F[cell].state, pv, ρ, X, Y, Sat, cell, N, has_other, phase_ix)
    end
end

function degrees_of_freedom_per_entity(model::SimulationModel{G,S}, v::TotalMasses) where {G<:Any,S<:MultiComponentSystem}
    number_of_components(model.system)
end

"""
Update total masses for two-phase compositional
"""
function two_phase_compositional_mass!(M, state, Φ, ρ, X, Y, S, cell, N, aqua::Val{false}, phase_ix)
    update_mass_two_phase_compositional!(M, state, Φ, ρ, X, Y, S, cell, phase_ix, N)
end

"""
Update total masses for two-phase compositional where another immiscible phase is present
"""
function two_phase_compositional_mass!(M, state, Φ, ρ, X, Y, S, cell, N, aqua::Val{true}, phase_ix)
    update_mass_two_phase_compositional!(M, state, Φ, ρ, X, Y, S, cell, phase_ix[2:end], N - 1)
    a, = phase_ix
    @inbounds M[end, cell] = ρ[a, cell]*S[a, cell]*Φ[cell]
end

function update_mass_two_phase_compositional!(M, state, Φ, ρ, X, Y, S, cell, phase_ix, N)
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
    @inbounds M_l = ρ[phase, cell] * S[phase, cell]
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
