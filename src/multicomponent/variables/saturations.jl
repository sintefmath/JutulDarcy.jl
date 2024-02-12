
@jutul_secondary function update_saturations!(Sat, s::Saturations, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ix) where {D, S<:MultiPhaseCompositionalSystemLV{E, T, O} where {E, T, O<:Nothing}}
    l, v = phase_indices(model.system)
    eos = model.system.equation_of_state
    @inbounds for i in ix
        S_l, S_v = phase_saturations(eos, Pressure[i], Temperature[i], FlashResults[i])
        Sat[l, i] = S_l
        Sat[v, i] = S_v
    end
end

@jutul_secondary function update_saturations!(Sat, s::Saturations, model::SimulationModel{D, S}, Pressure, Temperature, FlashResults, ImmiscibleSaturation, ix) where {D, S<:MultiPhaseCompositionalSystemLV}
    T = eltype(Sat)
    a, l, v = phase_indices(model.system)
    eos = model.system.equation_of_state
    @inbounds for i in ix
        S_other = ImmiscibleSaturation[i]
        fr = FlashResults[i]
        rem = one(T) - S_other + 0*MINIMUM_COMPOSITIONAL_SATURATION
        if fr.state == MultiComponentFlash.two_phase_lv
            S_v_pure = MultiComponentFlash.two_phase_vapor_saturation(eos, Pressure[i], Temperature[i], fr)
            S_v = rem*S_v_pure
        elseif fr.state == MultiComponentFlash.single_phase_v
            S_v = rem
        else
            S_v = zero(T)
        end
        S_l = rem - S_v
        Sat[l, i] = S_l
        Sat[v, i] = S_v
        Sat[a, i] = S_other
    end
end

struct SaturationsFromDensities <: PhaseVariables

end


@jutul_secondary function update_saturations!(Sat, s::SaturationsFromDensities, model::SimulationModel{D, S}, PhaseMassDensities, FlashResults, ix) where {D, S<:MultiPhaseCompositionalSystemLV{E, T, O} where {E, T, O<:Nothing}}
    l, v = phase_indices(model.system)
    @assert number_of_phases(model.system) == 2
    eos = model.system.equation_of_state
    props = map(x -> x.mw, eos.mixture.properties)
    @inbounds for i in ix
        flash = FlashResults[i]
        x = flash.liquid.mole_fractions
        y = flash.vapor.mole_fractions
        if flash.state == MultiComponentFlash.two_phase_lv
            # Calculate values that are proportional to volume and get
            # saturation from that.
            mass_liquid = 0.0
            mass_vapor = 0.0
            for c in eachindex(props)
                mw = props[i].mw
                mass_liquid += mw*x[i]
                mass_vapor += mw*y[i]
            end
            rho_l = PhaseMassDensities[l, i]
            rho_l = PhaseMassDensities[v, i]
            vol_liquid = mass_liquid/rho_l
            vol_vapor = mass_vapor/rho_v
            S_v = vol_vapor/(vol_liquid + vol_vapor)
        elseif flash.state == MultiComponentFlash.single_phase_v
            S_v = 1.0
        else
            S_v = 0.0
        end
        Sat[v, i] = S_v
        Sat[l, i] = 1.0 - S_v
    end
end
