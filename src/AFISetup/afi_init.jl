function setup_equilibrium_regions(d::AFIInputFile, model; regions = setup_region_map(d))
    model = reservoir_model(model)
    equils = find_records(d, "Equilibrium", "IX", steps = false, model = true, once = false)
    single_reg = length(equils) == 1
    return map(x -> setup_equilibrium_region(x.value, model, single_region = single_reg), equils)
end

function setup_equilibrium_region(equil, model; single_region = true)
    sys = model.system

    datum_depth = get(equil, "DatumDepth", 0.0)
    datum_pressure = get(equil, "DatumPressure", 0.0)
    woc_depth = get(equil, "WOCDepth", datum_depth)
    woc_pc = get(equil, "WOCPressure", 0.0)

    goc_depth = get(equil, "GOCDepth", datum_depth)
    goc_pc = get(equil, "GOCPressure", 0.0)

    if single_region
        cells = missing
    else
        region = equil["group"]
        error("Not implemented yet.")
    end
    rsvd = missing
    if JutulDarcy.has_disgas(sys)
        rs_tab = get(equil, "SolutionGORDepthTable", missing)
        if !ismissing(rs_tab)
            rsvd_tab = get_1d_interpolator(rs_tab["Depth"], rs_tab["SolutionGOR"])
            rsvd = z -> rsvd_tab(z)
        end
    end

    return EquilibriumRegion(model, datum_pressure;
        woc = woc_depth,
        goc = goc_depth,
        rs_vs_depth = rsvd,
        pc_woc = woc_pc,
        pc_goc = goc_pc,
        cells = cells
    )
end
