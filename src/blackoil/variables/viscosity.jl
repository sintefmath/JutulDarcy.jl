@jutul_secondary function update_as_secondary!(μ, ρ::DeckViscosity, model::StandardBlackOilModel, Pressure, Rs, Rv)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    nph, nc = size(μ)

    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    mb = minbatch(model.context, nc)
    @inbounds @batch minbatch = mb for i = 1:nc
        p = Pressure[i]
        rs = Rs[i]
        rv = Rv[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, rv, i)
    end
end

@jutul_secondary function update_as_secondary!(μ, ρ::DeckViscosity, model::VapoilBlackOilModel, Pressure, Rv)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    nph, nc = size(μ)

    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    mb = minbatch(model.context, nc)
    @inbounds @batch minbatch = mb for i = 1:nc
        p = Pressure[i]
        rv = Rv[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, i)
        μ[g, i] = viscosity(muG, reg, p, rv, i)
    end
end

@jutul_secondary function update_as_secondary!(μ, ρ::DeckViscosity, model::DisgasBlackOilModel, Pressure, Rs)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    nph, nc = size(μ)

    w, o, g = phase_indices(model.system)
    muW = pvt[w]
    muO = pvt[o]
    muG = pvt[g]

    mb = minbatch(model.context, nc)
    @inbounds @batch minbatch = mb for i = 1:nc
        p = Pressure[i]
        rs = Rs[i]
        μ[w, i] = viscosity(muW, reg, p, i)
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, i)
    end
end
