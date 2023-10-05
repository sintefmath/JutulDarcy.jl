
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
        rem = one(T) - S_other - MINIMUM_COMPOSITIONAL_SATURATION
        if fr.state == MultiComponentFlash.two_phase_lv
            S_v_pure = MultiComponentFlash.two_phase_vapor_saturation(eos, Pressure[i], Temperature[i], fr)
            S_v = rem*S_v_pure
        elseif fr.state == MultiComponentFlash.single_phase_v
            S_v = rem
        else
            S_v = zero(T)
        end
        # S_l, S_v = phase_saturations(eos, Pressure[i], Temperature[i], fr)
        S_l = one(T) - S_v - S_other
        Sat[l, i] = S_l
        Sat[v, i] = S_v
        Sat[a, i] = S_other
    end
end

