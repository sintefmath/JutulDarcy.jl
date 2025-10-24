function setup_saturation_variables(d::AFIInputFile, sys, reservoir; regions = setup_region_map(d))
    out = Dict{Symbol, Any}()
    out[:RelativePermeabilities] = setup_relperm(d, reservoir, sys, regions = regions)
    out[:CapillaryPressure] = setup_pc(d, reservoir, sys, regions = regions)
    return out
end

function setup_relperm(d, reservoir, sys; regions = setup_region_map(d))
    satfuns = find_records(d, "SaturationFunction", "IX", steps = false, model = true)

    phase_symb = phases_symbol(sys)
    present = phases_present(sys)
    krw = Dict()
    krg = Dict()
    krog = Dict()
    krow = Dict()
    drainage = get_region_value_and_map(reservoir, regions, "rock", "DRAINAGE_SATURATION_FUNCTION")
    imbibition = get_region_value_and_map(reservoir, regions, "rock", "IMBIBITION_SATURATION_FUNCTION")
    kr_regs = merge_saturation_regions(drainage, imbibition)

    for satfun in satfuns
        vals = satfun.value
        reg = satfun.value["region"]
        @warn reg
        if present.water
            krw[reg] = get_relperm(vals, :w, phase_symb, "WaterRelPermFunction")
            if present.oil
                krow[reg] = get_relperm(vals, :ow, phase_symb, "OilInWaterRelPermFunction")
            end
        end
        if present.gas
            krg[reg] = get_relperm(vals, :g, phase_symb, "GasRelPermFunction")
            if present.oil
                krog[reg] = get_relperm(vals, :og, phase_symb, "OilInGasRelPermFunction")
            end
        end
    end
    w = remap_to_tuple(krw, kr_regs)
    g = remap_to_tuple(krg, kr_regs)
    ow = remap_to_tuple(krow, kr_regs)
    og = remap_to_tuple(krog, kr_regs)

    hyst = find_records(d, "KilloughRelPermHysteresis", "IX", steps = false, model = true, once = true)
    if !isnothing(hyst)
        h = hyst.value
        # Drainage, imbibition, hysteresis?
        tol = get(h, "ModificationParameter", 0.1)
        nw_hyst = JutulDarcy.KilloughHysteresis(tol = tol)
        wetting_hyst = String(get(h, "WettingPhaseHysteresis", "DRAINAGE"))
        if wetting_hyst == "DRAINAGE"
            w_hyst = JutulDarcy.ImbibitionOnlyHysteresis()
        elseif wetting_hyst == "IMBIBITION"
            w_hyst = JutulDarcy.NoHysteresis()
        else
            wetting_hyst == "HYSTERESIS" || error("Unknown wetting phase hysteresis type: $wetting_hyst")
            w_hyst = nw_hyst
        end
        hysteresis_w = w_hyst
        hysteresis_ow = nw_hyst
        hysteresis_g = w_hyst
        hysteresis_og = nw_hyst
    else
        hysteresis_w = hysteresis_ow = hysteresis_g = hysteresis_og = JutulDarcy.NoHysteresis()
    end

    rock_mgr = find_records(d, "RockMgr", "IX", steps = false, model = true, once = true)
    hscale = get(rock_mgr.value, "HorizontalEndPointScaling", false)
    vscale = get(rock_mgr.value, "VerticalEndPointScaling", false)
    hscale == vscale || error("Only one of HorizontalEndPointScaling and VerticalEndPointScaling is set in RockMgr. Both must be true or both must be false in this current implementation.")
    if hscale || vscale
        if get(rock_mgr.value, "ThreePointScaling", false)
            scaling = JutulDarcy.ThreePointKrScale()
        else
            scaling = JutulDarcy.TwoPointKrScale()
        end
    else
        scaling = JutulDarcy.NoKrScale()
    end

    kr = JutulDarcy.ReservoirRelativePermeabilities(
        w = w,
        g = g,
        ow = ow,
        og = og,
        regions = kr_regs.value,
        hysteresis_w = hysteresis_w,
        hysteresis_ow = hysteresis_ow,
        hysteresis_g = hysteresis_g,
        hysteresis_og = hysteresis_og,
        scaling = scaling
    )
    return kr
end

function region_convert(region::Vector, reg_map = missing)
    has_reg_map = !ismissing(reg_map)

    minval = minimum(region)
    if has_reg_map
        minreg = minimum(values(reg_map))
        minval = min(minval, minreg)
    end
    offset = convert(Int, 1 - minval)
    if offset == 0 && eltype(region) == Int
        converted = region
    else
        converted = map(i -> Int(i + offset), region)
    end
    if offset != 0 && has_reg_map
        reg_map = copy(reg_map)
        for (k, v) in reg_map
            reg_map[k] = v + offset
        end
    end
    if has_reg_map
        out = (converted, reg_map)
    else
        out = converted
    end
    return out
end

function setup_pc(d, reservoir, sys; regions = setup_region_map(d))
    present = phases_present(sys)
    satfuns = find_records(d, "SaturationFunction", "IX", steps = false, model = true)
    pcow = Dict{String, Any}()
    pcog = Dict{String, Any}()

    has_ow = present.water && present.oil
    has_og = present.oil && present.gas

    function getpc(vals, name, sgn = 1.0)
        f = get_saturation_function(vals, "CapPressure", name)
        return get_1d_interpolator(f["Saturation"], sgn.*f["CapPressure"])
    end
    drainage = get_region_value_and_map(reservoir, regions, "rock", "DRAINAGE_SATURATION_FUNCTION")
    imbibition = get_region_value_and_map(reservoir, regions, "rock", "IMBIBITION_SATURATION_FUNCTION")
    regs = merge_saturation_regions(drainage, imbibition)
    for satfun in satfuns
        vals = satfun.value
        reg = satfun.value["region"]
        if has_ow
            pcow[reg] = getpc(vals, "OilWaterCapPressureFunction", -1)
        end
        if has_og
            pcog[reg] = getpc(vals, "GasOilCapPressureFunction")
        end
    end
    # TODO: Hysteresis, scaling
    pc = []
    if has_ow
        push!(pc, remap_to_tuple(pcow, regs))
    end
    if has_og
        push!(pc, remap_to_tuple(pcog, regs))
    end
    pc_f = JutulDarcy.SimpleCapillaryPressure(pc,
        regions = regs.value
    )
    return pc_f
end

function remap_to_tuple(d::Dict, regs)
    nkeys = length(keys(d))
    if nkeys == 0
        out = nothing
    else
        out = Vector{Any}(undef, nkeys)
        fill!(out, missing)
        for (k, v) in d
            out[regs.map[k]] = v
        end
        any(ismissing.(out)) && error("Not all regions mapped in saturation function remapping.")
        out = Tuple(out)
    end
    return out
end

function get_relperm(satfun, phase_label::Symbol, phases::Symbol, tab_label; throw = true)
    phase_label in (:w, :o, :ow, :og, :g) || error("Unsupported phase label for relperm: $phase_label")
    if haskey(satfun, "RelPerm")
        krtab = get_saturation_function(satfun, "RelPerm", tab_label, throw = throw)
        kr = krtab["RelPerm"]
        s = krtab["Saturation"]
    elseif haskey(satfun, "CoreyRelPerm")
        krtab = get_saturation_function(satfun, "CoreyRelPerm", tab_label, throw = throw)
        s, kr = brooks_corey_from_coefficients(krtab)
    end
    kr_obj = PhaseRelativePermeability(s, kr, label = phase_label)
    return kr_obj
end

function brooks_corey_from_coefficients(krtab)
    n = krtab["Exponent"]
    sr = get(krtab, "ResidualSaturation", 0.0)
    sconn = get(krtab, "ConnateSaturation", 0.0)
    krmax = get(krtab, "EndPointRelPerm", 1.0)
    s_max = get(krtab, "EndPointSaturation", 1.0 - sconn)
    s = collect(range(sconn, s_max, length = 50))
    kr = JutulDarcy.brooks_corey_relperm.(s, n = n, residual = sr, kr_max = krmax, residual_total = 1.0 - s_max)
    return (s, kr)
end

function get_saturation_function(satfun, type, label; throw = true)
    out = missing
    if haskey(satfun, label)
        table_name = satfun[label]
        return satfun[type][table_name]
    end
    if ismissing(out) && throw
        error("No $label for $type found in saturation function table.")
    end
    return out
end

function phases_present(sys)
    phases = JutulDarcy.get_phases(sys)
    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases
    return (water = has_water, oil = has_oil, gas = has_gas)
end

function phases_symbol(sys)
    present = phases_present(sys)
    has_water = present.water
    has_oil = present.oil
    has_gas = present.gas

    if has_water && has_oil && has_gas
        phases = :wog
    elseif has_water && has_oil
        phases = :wo
    elseif has_oil && has_gas
        phases = :og
    elseif has_water && has_gas
        phases = :wg
    elseif has_water
        phases = :w
    elseif has_oil
        phases = :o
    elseif has_gas
        phases = :g
    else
        error("No phases found in system")
    end
end

function get_region_value_and_map(reservoir, regions, regtype, regname)
    reg = get(regions[regtype], regname, missing)
    if ismissing(reg)
        out = missing
    else
        reg_no_label = only(keys(reg))
        reg_map = regions["family"][reg_no_label]
        reg_value = reservoir[Symbol(reg_no_label)]
        reg_value, reg_map = region_convert(reg_value, reg_map)
        out = (value = reg_value, map = reg_map)
    end
    return out
end

function merge_saturation_regions(drainage, imbibition)
    if ismissing(imbibition)
        merged = drainage
    else
        if !all(drainage.value .== imbibition.value)
            @warn "Drainage and imbibition saturation functions defined on different regions, cannot merge. Using drainage only."
        end
        new_map = copy(drainage.map)
        maxval = maximum(values(drainage.map))
        for (k, v) in imbibition.map
            if !haskey(new_map, k)
                new_map[k] = v + maxval
            end
        end
        merged = (value = drainage.value, map = merge(drainage.map, new_map))
    end
    return merged
end
