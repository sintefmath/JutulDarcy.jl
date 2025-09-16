function setup_equilibrium_regions(d::AFIInputFile, model; regions = setup_region_map(d))
    equils = find_records(d, "Equilibrium", "IX", steps = false, model = true, once = false)
    single_reg = length(equils) == 1
    function setup_equil(equil)
        datum_depth = get(equil, "DatumDepth", 0.0)
        datum_pressure = get(equil, "DatumPressure", 0.0)
        woc_depth = get(equil, "WOCDepth", datum_depth)
        woc_pc = get(equil, "WOCPressure", 0.0)

        goc_depth = get(equil, "GOCDepth", datum_depth)
        goc_pc = get(equil, "GOCPressure", 0.0)

        if single_reg
            cells = missing
        else
            region = equil["group"]
            error("Not implemented yet.")
        end
        return EquilibriumRegion(model, datum_pressure;
            woc = woc_depth,
            goc = goc_depth,
            pc_woc = woc_pc,
            pc_goc = goc_pc,
            cells = cells
        )
    end
    return map(x -> setup_equil(x.value), equils)
end
