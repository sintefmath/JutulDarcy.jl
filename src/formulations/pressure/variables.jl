





@jutul_secondary function update_pressure_factors(w, wdef::PressureReductionFactors, model, PhaseMassDensities, ix)
    for i in ix
        for ph in axes(w, 1)
            w[ph, i] = 1.0/PhaseMassDensities[ph, i]
        end
    end
end
