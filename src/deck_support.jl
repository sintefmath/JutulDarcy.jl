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


function Jutul.line_plot_data(model::SimulationModel, k::DeckPhaseVariables)
    deck_pvt_type(::DeckShrinkageFactors) = :shrinkage
    deck_pvt_type(::DeckDensity) = :density
    deck_pvt_type(::DeckViscosity) = :viscosity

    pvt_type = deck_pvt_type(k)
    has_reg = !isnothing(k.regions)
    s = collect(0:0.01:1)
    if has_reg
        nreg = length(k.relperms[1])
    else
        nreg = 1
    end
    nph = number_of_phases(model.system)
    phases = phase_names(model.system)
    data = Matrix{Any}(undef, nph, nreg)
    for r in 1:nreg
        for ph in 1:nph
            data[ph, r] = deck_function_plot_data(model, k.pvt[ph], ph, r, pvt_type)
        end
    end
    return data
end


function deck_function_plot_data(model, pvt, phase, reg, as_type)
    @info "" pvt as_type
    error("Not implemented for $(typeof(pvt)) as $as_type")
end

function deck_function_plot_data(model, pvt::Union{PVTW, PVDO, PVTW_EXTENDED, PVDG}, phase, reg, as_type)
    p = collect(range(1e5, 1000e5, 100))
    phase_name = phase_names(model.system)[phase]
    if as_type == :shrinkage
        F = (pressure) -> shrinkage(pvt.tab[reg], pressure)
        title = "Phase Shrinkage"
        yl = "1/B"
    elseif as_type == :viscosity
        F = (pressure) -> viscosity(pvt.tab[reg], pressure)/1e-3
        title = "Phase Viscosity"
        yl = "μ [cP]"
    else
        @assert as_type == :density
        rhoS = reference_densities(model.system)[phase]
        F = (pressure) -> rhoS*shrinkage(pvt.tab[reg], pressure)
        title = "Phase Mass Density"
        yl = "ρ [kg/m^3]"
    end
    return Jutul.JutulLinePlotData(p./1e5, F.(p), labels = phase_name, title = title, xlabel = "Pressure [bar]", ylabel = yl)
end
