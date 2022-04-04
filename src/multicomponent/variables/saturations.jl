
@jutul_secondary function update_as_secondary!(Sat, s::Saturations, model::SimulationModel{D, S}, param, Pressure, Temperature, FlashResults) where {D, S<:MultiPhaseCompositionalSystemLV{E, T, O} where {E, T, O<:Nothing}}
    tb = thread_batch(model.context)
    l, v = phase_indices(model.system)
    eos = model.system.equation_of_state
    @inbounds @batch minbatch = tb for i in 1:size(Sat, 2)
        S_l, S_v = phase_saturations(eos, Pressure[i], Temperature[i], FlashResults[i])
        Sat[l, i] = S_l
        Sat[v, i] = S_v
    end
end

@jutul_secondary function update_as_secondary!(Sat, s::Saturations, model::SimulationModel{D, S}, param, Pressure, Temperature, FlashResults, ImmiscibleSaturation) where {D, S<:MultiPhaseCompositionalSystemLV}
    T = eltype(Sat)
    a, l, v = phase_indices(model.system)
    tb = thread_batch(model.context)
    eos = model.system.equation_of_state
    @inbounds @batch minbatch = tb for i in 1:size(Sat, 2)
        S_other = ImmiscibleSaturation[i]
        rem = one(T) - S_other
        S_l, S_v = phase_saturations(eos, Pressure[i], Temperature[i], FlashResults[i])
        Sat[l, i] = rem*S_l
        Sat[v, i] = rem*S_v
        Sat[a, i] = S_other
    end
end

