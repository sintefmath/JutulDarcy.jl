
@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::DisgasBlackOilModel, Pressure, Rs)
    pvt, reg = ρ.pvt, ρ.regions
    nph, nc = size(b)
    tb = minbatch(model.context, nc)

    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds @batch minbatch = tb for i in 1:nc
        p = Pressure[i]
        rs = Rs[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, i)
    end
end

# Shrinkage factors for all three cases
@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::StandardBlackOilModel, Pressure, Rs, Rv)
    pvt, reg = ρ.pvt, ρ.regions
    nph, nc = size(b)
    tb = minbatch(model.context, nc)

    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds @batch minbatch = tb for i in 1:nc
        p = Pressure[i]
        rv = Rv[i]
        rs = Rs[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, rv, i)
    end
end

@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model::VapoilBlackOilModel, Pressure, Rv)
    pvt, reg = ρ.pvt, ρ.regions
    nph, nc = size(b)
    tb = minbatch(model.context, nc)

    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds @batch minbatch = tb for i in 1:nc
        p = Pressure[i]
        rv = Rv[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, i)
        b[g, i] = shrinkage(bG, reg, p, rv, i)
    end
end
