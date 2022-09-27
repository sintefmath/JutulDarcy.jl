export DeckViscosity, DeckShrinkage

@jutul_secondary function update_as_secondary!(mu, μ::DeckViscosity, model, Pressure)
    pvt, reg = μ.pvt, μ.regions
    if false
        @tullio mu[ph, i] = viscosity(pvt[ph], reg, Pressure[i], i)
    else
        tb = minbatch(model.context)
        nph, nc = size(mu)
        for ph in 1:nph
            pvt_ph = pvt[ph]
            @batch minbatch = tb for i in 1:nc
                p = Pressure[i]
                @inbounds mu[ph, i] = viscosity(pvt_ph, reg, p, i)
            end
        end
    end
end

@jutul_secondary function update_as_secondary!(rho, ρ::DeckDensity, model, Pressure)
    rhos = reference_densities(model.system)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    tb = minbatch(model.context)
    nph, nc = size(rho)
    for ph in 1:nph
        rhos_ph = rhos[ph]
        pvt_ph = pvt[ph]
        @batch minbatch = tb for i in 1:nc
            p = Pressure[i]
            @inbounds rho[ph, i] = rhos_ph*shrinkage(pvt_ph, reg, p, i)
        end
    end
end


@jutul_secondary function update_as_secondary!(b, ρ::DeckShrinkageFactors, model, Pressure)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    tb = minbatch(model.context)
    nph, nc = size(b)
    for ph in 1:nph
        pvt_ph = pvt[ph]
        @batch minbatch = tb for i in 1:nc
            p = Pressure[i]
            @inbounds b[ph, i] = shrinkage(pvt_ph, reg, p, i)
        end
    end
end

@jutul_secondary function update_as_secondary!(pv, Φ::LinearlyCompressiblePoreVolume, model, Pressure, StaticFluidVolume)
    c_r = Φ.expansion
    p_r = Φ.reference_pressure
    @inbounds for i in eachindex(pv)
        p = Pressure[i]
        x = c_r*(p-p_r)
        mult = 1.0 + x + 0.5*(x^2)
        pv[i] = StaticFluidVolume[i]*mult
    end
end
