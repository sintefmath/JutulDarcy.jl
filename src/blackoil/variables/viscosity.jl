@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckViscosity, model::StandardBlackOilModel, Pressure, Rs, Rv, ix)
    pvt, reg = ρ.pvt, ρ.regions
    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    @inbounds for i in ix
        p = Pressure[i]
        rs = Rs[i]
        rv = Rv[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, rv, i)
    end
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckViscosity, model::VapoilBlackOilModel, Pressure, Rv, ix)
    pvt, reg = ρ.pvt, ρ.regions
    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    @inbounds for i in ix
        p = Pressure[i]
        rv = Rv[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, i)
        μ[g, i] = viscosity(muG, reg, p, rv, i)
    end
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckViscosity, model::DisgasBlackOilModel, Pressure, Rs, ix)
    pvt, reg = ρ.pvt, ρ.regions
    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]
    @inbounds for i in ix
        p = Pressure[i]
        rs = Rs[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, i)
    end
end
