function setup_equilibrium_regions(d::AFIInputFile, model; regions = setup_region_map(d))
    model = reservoir_model(model)
    equils = find_records(d, "Equilibrium", "IX", steps = false, model = true, once = false)
    equil = merge_records(equils).value
    return map(
        reg -> setup_equilibrium_region(equil[reg], model),
        collect(keys(equil))
    )
end

function merge_records(vals)
    I = firstindex(vals)
    kw = vals[I].keyword
    a = vals[I].value
    for i in eachindex(vals)
        if i == I
            continue
        end
        @assert vals[i].keyword == kw "Cannot merge dissimilar keywords: '$kw' and '$(vals[i].keyword)'"
        b = vals[i].value
        a = merge_records(a, b)
    end
    return (value = a, keyword = kw)
end

function merge_records(a, b)
    keysa = keys(a)
    keysb = keys(b)
    out = Dict{String, Any}()
    for k in intersect(keysa, keysb)
        out[k] = merge(a[k], b[k])
    end
    for ka in setdiff(keysa, keysb)
        out[ka] = a[ka]
    end
    for kb in setdiff(keysb, keysa)
        out[kb] = b[kb]
    end
    return out
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
    if JutulDarcy.model_is_thermal(model)
        tdt = equil["TemperatureDepthTable"]
        T_d = tdt["Temperature"] .+ 273.15
        I_T = get_1d_interpolator(
            tdt["Depth"],
            T_d
        )
        temperature_vs_depth = z -> I_T(z)
    else
        temperature_vs_depth = missing
    end

    return EquilibriumRegion(model, datum_pressure;
        woc = woc_depth,
        goc = goc_depth,
        rs_vs_depth = rsvd,
        pc_woc = woc_pc,
        pc_goc = goc_pc,
        temperature_vs_depth = temperature_vs_depth,
        cells = cells
    )
end
