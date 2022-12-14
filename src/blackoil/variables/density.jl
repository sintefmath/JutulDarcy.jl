# Density for all three cases
@jutul_secondary function update_deck_density!(rho, m::DeckDensity, model::StandardBlackOilModel, Rs, Rv, ShrinkageFactors, ix)
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

@jutul_secondary function update_deck_density!(rho, m::DeckDensity, model::VapoilBlackOilModel, Rv, ShrinkageFactors, ix)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    for i in ix
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*rhoOS
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_deck_density!(rho, m::DeckDensity, model::DisgasBlackOilModel, Rs, ShrinkageFactors, ix)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    for i in ix
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*rhoGS
    end
end
