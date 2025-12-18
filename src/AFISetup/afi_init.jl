function setup_equilibrium_regions(d::AFIInputFile, model; regions = setup_region_map(d))
    model = reservoir_model(model)
    reservoir = reservoir_domain(model)
    equils = find_records(d, "Equilibrium", "IX", steps = false, model = true, once = false)
    regmap = get_region_value_and_map(reservoir, regions, "equilibrium", "EQUILIBRIUM_EQL_MODEL")
    regs_eql = regions["equilibrium"]["EQUILIBRIUM_EQL_MODEL"]
    equil = merge_records(equils).value
    return map(
        reg -> setup_equilibrium_region(equil[reg], model, cells = equil_model_region_to_cells(regmap, regs_eql, reg)),
        collect(keys(equil))
    )
end

function equil_model_region_to_cells(regmap, regs_eql, model_reg)
    equil_reg = missing
    for (k, v) in pairs(regs_eql)
        for el in v
            if el.model == model_reg
                equil_reg = el.region
                break
            end
        end
    end
    !ismissing(equil_reg) || error("No equilibrium region found for model region $model_reg")
    return findall(isequal(regmap.map[equil_reg]), regmap.value)
end

function setup_equilibrium_region(equil, model; cells = missing)
    sys = model.system

    datum_depth = get(equil, "DatumDepth", 0.0)
    datum_pressure = get(equil, "DatumPressure", 0.0)
    woc_depth = get(equil, "WOCDepth", datum_depth)
    woc_pc = get(equil, "WOCPressure", 0.0)

    goc_depth = get(equil, "GOCDepth", datum_depth)
    goc_pc = get(equil, "GOCPressure", 0.0)

    rsvd = missing
    if JutulDarcy.has_disgas(sys)
        rs_tab = get(equil, "SolutionGORDepthTable", missing)
        if !ismissing(rs_tab)
            rsvd_tab = get_1d_interpolator(rs_tab["Depth"], rs_tab["SolutionGOR"])
            rsvd = z -> rsvd_tab(z)
        end
    end
    if JutulDarcy.model_is_thermal(model)
        tdt = equil["TemperatureDepthTable"]
        T_d = tdt["Temperature"] .+ 273.15
        I_T = get_1d_interpolator(
            tdt["Depth"],
            T_d,
            cap_endpoints = false
        )
        temperature_vs_depth = z -> I_T(z)
    else
        temperature_vs_depth = missing
    end

    return EquilibriumRegion(model, datum_pressure, datum_depth;
        woc = woc_depth,
        goc = goc_depth,
        rs_vs_depth = rsvd,
        pc_woc = woc_pc,
        pc_goc = goc_pc,
        temperature_vs_depth = temperature_vs_depth,
        cells = cells
    )
end
