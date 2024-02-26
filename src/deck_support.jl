@jutul_secondary function update_deck_viscosity!(mu, μ::DeckPhaseViscosities, model, Pressure, ix)
    pvt, reg = μ.pvt, μ.regions
    @inbounds for i in ix
        p = Pressure[i]
        for ph in axes(mu, 1)
            pvt_ph = pvt[ph]
            @inbounds mu[ph, i] = viscosity(pvt_ph, reg, p, i)
        end
    end
end

@jutul_secondary function update_deck_viscosity!(mu, μ::DeckPhaseViscosities{<:Any, Ttab, <:Any}, model, Pressure, Temperature, ix) where Ttab<:DeckThermalViscosityTable
    pvt, reg = μ.pvt, μ.regions
    for i in ix
        r_i = region(μ, i)
        p = Pressure[i]
        T = Temperature[i]
        for ph in axes(mu, 1)
            pvt_ph = pvt[ph]
            pvt_thermal = table_by_region(μ.thermal.visc_tab[ph], r_i)
            p_ref = table_by_region(μ.thermal.p_ref[ph], r_i)

            mu_p = viscosity(pvt_ph, reg, p, i)
            mu_ref = viscosity(pvt_ph, reg, p_ref, i)
            mu_thermal = pvt_thermal(T)
            mu[ph, i] = mu_thermal*(mu_p/mu_ref)
        end
    end
end

@jutul_secondary function update_deck_density!(rho, ρ::DeckPhaseMassDensities, model, Pressure, ix)
    rhos = reference_densities(model.system)
    pvt, reg = ρ.pvt, ρ.regions
    # Note immiscible assumption
    @inbounds for ph in axes(rho, 1)
        rhos_ph = rhos[ph]
        pvt_ph = pvt[ph]
        @inbounds for i in ix
            p = Pressure[i]
            rho[ph, i] = rhos_ph*shrinkage(pvt_ph, reg, p, i)
        end
    end
end

@jutul_secondary function update_deck_density!(rho, ρ::DeckPhaseMassDensities{<:Any, <:WATDENT, <:Any}, model, Pressure, Temperature, ix)
    rhos = reference_densities(model.system)
    pvt, reg = ρ.pvt, ρ.regions
    phases = get_phases(model.system)
    # Note immiscible assumption
    for i in ix
        r_i = region(ρ, i)
        p = Pressure[i]
        T = Temperature[i]
        for ph in axes(rho, 1)
            rhos_ph = rhos[ph]
            pvt_ph = pvt[ph]
            if phases[ph] == AqueousPhase()
                T_ref, c1, c2 = ρ.watdent.tab[r_i]
                pvtw = pvt_ph.tab[r_i]
                p_ref = pvtw.p_ref
                B_pref = 1.0/shrinkage(pvt_ph, reg, p_ref, i)
                Δp = pvtw.b_c*(p - p_ref)
                ΔT = T - T_ref
                B_w = B_pref*(1.0 - Δp)*(1.0 + c1*ΔT + c2*ΔT^2)
                rho[ph, i] = rhos_ph/B_w
            else
                rho[ph, i] = rhos_ph*shrinkage(pvt_ph, reg, p, i)
            end
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

function Jutul.line_plot_data(model::SimulationModel, k::LinearlyCompressiblePoreVolume)
    p = collect(range(1e5, 1000e5, 1000))
    v = ones(size(p))
    y = similar(v)
    update_pore_volume!(y, k, model, p, v, eachindex(y))
    Jutul.JutulLinePlotData(p./1e5, y, labels = "Φ", title = "Rock pore volume expansion", xlabel = "Pressure [bar]", ylabel = "Pore-volume multiplier")
end

function Jutul.line_plot_data(model::SimulationModel, k::DeckPhaseVariables)
    deck_pvt_type(::DeckShrinkageFactors) = :shrinkage
    deck_pvt_type(::DeckPhaseMassDensities) = :density
    deck_pvt_type(::DeckPhaseViscosities) = :viscosity

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
            x = deck_function_plot_data(model, k.pvt[ph], ph, r, pvt_type)
            if ismissing(x)
                continue
            end
            data[ph, r] = x
        end
    end
    return data
end


function deck_function_plot_data(model, pvt, phase, reg, as_type)
    @warn "Plotting not implemented for $(typeof(pvt)) as $as_type"
    return missing
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

function deck_function_plot_data(model, pvt::Union{PVTG, PVTO}, phase, reg, as_type)
    sys = model.system
    phase_name = phase_names(sys)[phase]
    rhoS = reference_densities(sys)
    rhoS_self = rhoS[phase]
    
    ix = phase_indices(sys)
    if length(ix) == 3
        a, l, v = ix
    else
        l, v = ix
    end
    if pvt isa PVTG
        sat_fn = model.system.rv_max
        @assert phase == v
        rhoS_other = rhoS[l]
    else
        sat_fn = model.system.rs_max
        @assert phase == l
        rhoS_other = rhoS[v]
    end
    if as_type == :shrinkage
        F = (pressure, r) -> shrinkage(pvt, reg, pressure, r, 1)
        title = "Phase Shrinkage"
        yl = "1/B"
    elseif as_type == :viscosity
        F = (pressure, r) -> viscosity(pvt, reg, pressure, r, 1)/1e-3
        title = "Phase Viscosity"
        yl = "μ [cP]"
    else
        @assert as_type == :density
        F = (pressure, r) -> (rhoS_self + r*rhoS_other)*shrinkage(pvt, reg, pressure, r, 1)
        title = "Phase Mass Density"
        yl = "ρ [kg/m^3]"
    end
    np = 100
    nsat = 10
    data = Vector{Float64}()
    p = Vector{Float64}()

    sizehint!(data, np*nsat)
    sizehint!(p, np*nsat)

    p_max = maximum(sat_fn.X)
    for r_ix in 1:nsat
        for p_i in range(1e5, p_max, np)
            r = (r_ix-1)*sat_fn(p_i)/(nsat-1)
            push!(p, p_i)
            push!(data, F(p_i, r))
        end
        push!(p, NaN)
        push!(data, NaN)
    end
    return Jutul.JutulLinePlotData(p./1e5, data, labels = phase_name, title = title, xlabel = "Pressure [bar]", ylabel = yl)
end
