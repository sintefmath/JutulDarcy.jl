
@jutul_secondary function update_deck_shrinkage!(b, ρ::DeckShrinkageFactors, model::DisgasBlackOilModel,
                                             Pressure, Rs, ix)
    pvt, reg = ρ.pvt, ρ.regions
    sys = model.system
    has_water = has_other_phase(sys)
    if has_water
        w, o, g = phase_indices(sys)
        bW = pvt[w]
    else
        o, g = phase_indices(sys)
    end
    bO = pvt[o]
    bG = pvt[g]
    @inbounds for i in ix
        p = Pressure[i]
        rs = Rs[i]
        if has_water
            b[w, i] = shrinkage(bW, reg, p, i)
        end
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, i)
    end
end

# Shrinkage factors for all three cases
@jutul_secondary function update_deck_shrinkage!(b, ρ::DeckShrinkageFactors, model::StandardBlackOilModel, Pressure, Rs, Rv, ix)
    pvt, reg = ρ.pvt, ρ.regions
    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds for i in ix
        p = Pressure[i]
        rv = Rv[i]
        rs = Rs[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, rs, i)
        b[g, i] = shrinkage(bG, reg, p, rv, i)
    end
end

@jutul_secondary function update_deck_shrinkage!(b, ρ::DeckShrinkageFactors, model::VapoilBlackOilModel, Pressure, Rv, ix)
    pvt, reg = ρ.pvt, ρ.regions
    w, o, g = phase_indices(model.system)
    bO = pvt[o]
    bG = pvt[g]
    bW = pvt[w]
    @inbounds for i in ix
        p = Pressure[i]
        rv = Rv[i]
        b[w, i] = shrinkage(bW, reg, p, i)
        b[o, i] = shrinkage(bO, reg, p, i)
        b[g, i] = shrinkage(bG, reg, p, rv, i)
    end
end
