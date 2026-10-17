function update_blackoil_viscosities!(μ, ρ::DeckPhaseViscosities, model::StandardBlackOilModel, Pressure, Rs, Rv, ix)
    pvt, reg = ρ.pvt, ρ.regions
    has_wat = has_other_phase(model.system)
    if has_wat
        w, o, g = phase_indices(model.system)
        muW = pvt[w]
    else
        o, g = phase_indices(model.system)
    end
    muO = pvt[o]
    muG = pvt[g]

    @inbounds for i in ix
        p = Pressure[i]
        rs = Rs[i]
        rv = Rv[i]
        if has_wat
            μ[w, i] = viscosity(muW, reg, p, i)
        end
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, rv, i)
    end
end

function update_blackoil_viscosities!(μ, ρ::DeckPhaseViscosities, model::VapoilBlackOilModel, Pressure, Rv, ix)
    pvt, reg = ρ.pvt, ρ.regions
    has_wat = has_other_phase(model.system)
    if has_wat
        w, o, g = phase_indices(model.system)
        muW = pvt[w]
    else
        o, g = phase_indices(model.system)
    end
    muO = pvt[o]
    muG = pvt[g]

    @inbounds for i in ix
        p = Pressure[i]
        rv = Rv[i]
        if has_wat
            μ[w, i] = viscosity(muW, reg, p, i)
        end
        μ[o, i] = viscosity(muO, reg, p, i)
        μ[g, i] = viscosity(muG, reg, p, rv, i)
    end
end

function update_blackoil_viscosities!(μ, ρ::DeckPhaseViscosities, model::DisgasBlackOilModel, Pressure, Rs, ix)
    pvt, reg = ρ.pvt, ρ.regions
    has_wat = has_other_phase(model.system)
    if has_wat
        w, o, g = phase_indices(model.system)
        muW = pvt[w]
    else
        o, g = phase_indices(model.system)
    end
    muO = pvt[o]
    muG = pvt[g]
    @inbounds for i in ix
        p = Pressure[i]
        rs = Rs[i]
        if has_wat
            μ[w, i] = viscosity(muW, reg, p, i)
        end
        μ[o, i] = viscosity(muO, reg, p, rs, i)
        μ[g, i] = viscosity(muG, reg, p, i)
    end
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckPhaseViscosities, model::StandardBlackOilModel, Pressure, Rs, Rv, ix)
    update_blackoil_viscosities!(μ, ρ, model, Pressure, Rs, Rv, ix)
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckPhaseViscosities, model::VapoilBlackOilModel, Pressure, Rv, ix)
    update_blackoil_viscosities!(μ, ρ, model, Pressure, Rv, ix)
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckPhaseViscosities, model::DisgasBlackOilModel, Pressure, Rs, ix)
    update_blackoil_viscosities!(μ, ρ, model, Pressure, Rs, ix)
end

# Resolve both the thermal-table and black-oil specializations explicitly so
# that dependencies include temperature as well as the solution ratios.
@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckPhaseViscosities{<:Any, <:DeckThermalViscosityTable}, model::DisgasBlackOilModel, Pressure, Temperature, Rs, ix)
    update_blackoil_viscosities!(μ, ρ, model, Pressure, Rs, ix)
    apply_blackoil_thermal_viscosity!(μ, ρ, model, Temperature, Rs, nothing, ix)
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckPhaseViscosities{<:Any, <:DeckThermalViscosityTable}, model::VapoilBlackOilModel, Pressure, Temperature, Rv, ix)
    update_blackoil_viscosities!(μ, ρ, model, Pressure, Rv, ix)
    apply_blackoil_thermal_viscosity!(μ, ρ, model, Temperature, nothing, Rv, ix)
end

@jutul_secondary function update_deck_viscosity!(μ, ρ::DeckPhaseViscosities{<:Any, <:DeckThermalViscosityTable}, model::StandardBlackOilModel, Pressure, Temperature, Rs, Rv, ix)
    update_blackoil_viscosities!(μ, ρ, model, Pressure, Rs, Rv, ix)
    apply_blackoil_thermal_viscosity!(μ, ρ, model, Temperature, Rs, Rv, ix)
end

function apply_blackoil_thermal_viscosity!(μ, variable, model, Temperature, Rs, Rv, ix)
    foreach_phase_pvt(variable.pvt) do phase_id, pvt
        ph = phase_index(phase_id)
        phase = phase_tuple_entry(get_phases(model.system), phase_id)
        tables = phase_tuple_entry(variable.thermal.visc_tab, phase_id)
        if isnothing(tables)
            return nothing
        end
        for i in ix
            r = region(variable, i)
            thermal_mu = evaluate_table_by_region(tables, r, Temperature[i])
            p_ref = table_by_region(variable.thermal.p_ref, r)
            if isfinite(p_ref)
                reg = variable.regions
                if phase == LiquidPhase() && !isnothing(Rs)
                    rs_ref = table_by_region(variable.thermal.rs_ref, r)
                    if isfinite(rs_ref)
                        rs = rs_ref
                    else
                        rs = Rs[i]
                    end
                    mu_ref = viscosity(pvt, reg, p_ref, rs, i)
                elseif phase == VaporPhase() && !isnothing(Rv)
                    mu_ref = viscosity(pvt, reg, p_ref, Rv[i], i)
                else
                    mu_ref = viscosity(pvt, reg, p_ref, i)
                end
                thermal_mu *= μ[ph, i]/mu_ref
            end
            μ[ph, i] = thermal_mu
        end
    end
    return μ
end
