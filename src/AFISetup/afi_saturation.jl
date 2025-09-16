function setup_saturation_variables(d::AFIInputFile, sys, reservoir; regions = setup_region_map(d))
    out = Dict{Symbol, Any}()
    out[:RelativePermeabilities] = setup_relperm(d, reservoir, sys, regions = regions)
    out[:CapillaryPressure] = setup_pc(d, reservoir, sys, regions = regions)
    return out
end

function setup_relperm(d, reservoir, sys; regions = setup_region_map(d))
    phase_symb = phases_symbol(sys)
    present = phases_present(sys)
    satfuns = find_records(d, "SaturationFunction", "IX", steps = false, model = true)
    krw = Dict()
    krg = Dict()
    krog = Dict()
    krow = Dict()
    for satfun in map(x -> x, satfuns)
        vals = satfun.value
        reg = satfun.value["region"]
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
    drainage_satnum, drain_reg_map = get_region_value_and_map(reservoir, regions, "rock", "DRAINAGE_SATURATION_FUNCTION")
    # TODO: Hysteresis, scaling
    kr = JutulDarcy.ReservoirRelativePermeabilities(
        w = remap_to_tuple(krw, drain_reg_map),
        g = remap_to_tuple(krg, drain_reg_map),
        ow = remap_to_tuple(krow, drain_reg_map),
        og = remap_to_tuple(krog, drain_reg_map),
        regions = drainage_satnum
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
    for satfun in map(x -> x, satfuns)
        vals = satfun.value
        reg = satfun.value["region"]
        if has_ow
            pcow[reg] = get_saturation_function(vals, "CapPressure", "OilWaterCapPressureFunction")
        end
        if has_og
            pcog[reg] = get_saturation_function(vals, "CapPressure", "GasOilCapPressureFunction")
        end
    end
    drainage_satnum, drain_reg_map = get_region_value_and_map(reservoir, regions, "rock", "DRAINAGE_SATURATION_FUNCTION")
    # TODO: Hysteresis, scaling
    pc = []
    if has_ow
        push!(pc, remap_to_tuple(pcow, drain_reg_map))
    end
    if has_og
        push!(pc, remap_to_tuple(pcog, drain_reg_map))
    end
    pc_f = JutulDarcy.SimpleCapillaryPressure(pc,
        regions = drainage_satnum
    )
    return pc_f
end

function remap_to_tuple(d::Dict, mapper)
    nkeys = length(keys(d))
    if nkeys == 0
        out = nothing
    else
        out = Vector{Any}(undef, nkeys)
        fill!(out, missing)
        for (k, v) in d
            out[mapper[k]] = v
        end
        any(ismissing(out)) && error("Not all regions mapped in saturation function remapping.")
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
        n = krtab["Exponent"]
        sr = get(krtab, "ResidualSaturation", 0.0)
        sconn = get(krtab, "ConnateSaturation", 0.0)
        krmax = get(krtab, "EndPointRelPerm", 1.0)
        s_max = get(krtab, "EndPointSaturation", 1.0 - sconn)
        s = collect(range(sconn, s_max, length = 50))
        kr = JutulDarcy.brooks_corey_relperm.(s, n = n, residual = sr, kr_max = krmax, residual_total = 1.0 - s_max)
    end
    return PhaseRelativePermeability(s, kr, label = phase_label)
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
    reg = regions[regtype][regname]
    reg_no_label = only(keys(reg))
    # drain_reg_no = only(values(drain_reg))

    reg_map = regions["family"][reg_no_label]
    reg_value = reservoir[Symbol(reg_no_label)]
    reg_value, reg_map = region_convert(reg_value, reg_map)
    return (reg_value, reg_map)
end
