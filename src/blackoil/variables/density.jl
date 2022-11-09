# Density for all three cases
@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::StandardBlackOilModel, Rs, Rv, ShrinkageFactors)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    n = size(rho, 2)
    mb = minbatch(model.context, n)
    @inbounds @batch minbatch = mb for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::VapoilBlackOilModel, Rv, ShrinkageFactors)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    n = size(rho, 2)
    mb = minbatch(model.context, n)
    @inbounds @batch minbatch = mb for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*rhoOS
        rho[g, i] = b[g, i]*(rhoGS + Rv[i]*rhoOS)
    end
end

@jutul_secondary function update_as_secondary!(rho, m::DeckDensity, model::DisgasBlackOilModel, Rs, ShrinkageFactors)
    b = ShrinkageFactors
    sys = model.system
    w, o, g = phase_indices(sys)
    rhoWS, rhoOS, rhoGS = reference_densities(sys)
    n = size(rho, 2)
    mb = minbatch(model.context, n)
    @inbounds @batch minbatch = mb for i = 1:n
        rho[w, i] = b[w, i]*rhoWS
        rho[o, i] = b[o, i]*(rhoOS + Rs[i]*rhoGS)
        rho[g, i] = b[g, i]*rhoGS
    end
end
