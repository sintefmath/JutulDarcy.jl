# Density for all three cases
@jutul_secondary function update_deck_density!(rho, m::DeckPhaseMassDensities, model::StandardBlackOilModel, Rs, Rv, ShrinkageFactors, ix)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    @inbounds for i in ix
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_deck_density!(rho, m::DeckPhaseMassDensities, model::VapoilBlackOilModel, Rv, ShrinkageFactors, ix)
    b = ShrinkageFactors
    sys = model.system
    has_water = has_other_phase(sys)
    if has_water
        w, o, g = phase_indices(sys)
        rhoWS, rhoOS, rhoGS = reference_densities(sys)
    else
        o, g = phase_indices(sys)
        rhoOS, rhoGS = reference_densities(sys)
    end
    for i in ix
        if has_water
            rho[w, i] = b[w, i]*rhoWS
        end
        rho[o, i] = b[o, i]*rhoOS
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_deck_density!(rho, m::DeckPhaseMassDensities, model::DisgasBlackOilModel, Rs, ShrinkageFactors, ix)
    b = ShrinkageFactors
    sys = model.system
    has_water = has_other_phase(sys)
    if has_water
        w, o, g = phase_indices(sys)
        rhoWS, rhoOS, rhoGS = reference_densities(sys)
    else
        o, g = phase_indices(sys)
        rhoOS, rhoGS = reference_densities(sys)
    end
    for i in ix
        if has_water
            rho[w, i] = b[w, i]*rhoWS
        end
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*rhoGS
    end
end
