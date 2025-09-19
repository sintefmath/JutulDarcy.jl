function setup_equilibrium_regions(d::AFIInputFile, reservoir, sys; regions = setup_region_map(d))
    equils = find_records(d, "Equilibrium", "IX", steps = false, model = true, once = false)
    for equil in equils
        equil = equil.value
        @info "??" equil
        region = equil["group"]
        datum_depth = get(equil, "DatumDepth", 0.0)
        datum_pressure = get(equil, "DatumPressure", 0.0)
        woc_depth = get(equil, "WOCDepth", datum_depth)
        woc_pc = get(equil, "WOCPressure", 0.0)

        goc_depth = get(equil, "GOCDepth", datum_depth)
        goc_pc = get(equil, "GOCPressure", 0.0)

        
    end
end
