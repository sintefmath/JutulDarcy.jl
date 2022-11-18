export DeckViscosity, DeckShrinkage

@jutul_secondary function update_deck_viscosity!(mu, μ::DeckViscosity, model, Pressure, ix)
    pvt, reg = μ.pvt, μ.regions
    nph = size(mu, 1)
    for ph in 1:nph
        pvt_ph = pvt[ph]
        for i in ix
            p = Pressure[i]
            @inbounds mu[ph, i] = viscosity(pvt_ph, reg, p, i)
        end
    end
end

@jutul_secondary function update_deck_density!(rho, ρ::DeckDensity, model, Pressure, ix)
    rhos = reference_densities(model.system)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    nph, nc = size(rho)
    for ph in 1:nph
        rhos_ph = rhos[ph]
        pvt_ph = pvt[ph]
        for i in ix
            p = Pressure[i]
            @inbounds rho[ph, i] = rhos_ph*shrinkage(pvt_ph, reg, p, i)
        end
    end
end


@jutul_secondary function update_deck_shrinkage!(b, ρ::DeckShrinkageFactors, model, Pressure, ix)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    nph, = size(b, 1)
    for ph in 1:nph
        pvt_ph = pvt[ph]
        for i in ix
            p = Pressure[i]
            @inbounds b[ph, i] = shrinkage(pvt_ph, reg, p, i)
        end
    end
end

@jutul_secondary function update_pore_volume!(pv, Φ::LinearlyCompressiblePoreVolume, model, Pressure, StaticFluidVolume, ix)
    c_r = Φ.expansion
    p_r = Φ.reference_pressure
    @inbounds for i in ix
        p = Pressure[i]
        x = c_r*(p-p_r)
        mult = 1.0 + x + 0.5*(x^2)
        pv[i] = StaticFluidVolume[i]*mult
    end
end
