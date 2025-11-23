function setup_saturation_variables(d::AFIInputFile, sys, reservoir;
        regions = setup_region_map(d),
        disable_endscale = false,
        disable_hysteresis = false,
    )
    out = Dict{Symbol, Any}()
    if length(JutulDarcy.get_phases(sys)) > 1
        out[:RelativePermeabilities] = setup_relperm(d, reservoir, sys,
            regions = regions,
            disable_endscale = disable_endscale,
            disable_hysteresis = disable_hysteresis
        )
        out[:CapillaryPressure] = setup_pc(d, reservoir, sys, regions = regions)
    end
    return out
end

function setup_relperm(d, reservoir, sys;
        regions = setup_region_map(d),
        disable_endscale = false,
        disable_hysteresis = false
    )
    satfuns = find_records(d, "SaturationFunction", "IX", steps = false, model = true)

    phase_symb = phases_symbol(sys)
    present = phases_present(sys)
    length(present) > 1 || error("Relative permeability setup called for single-phase system.")
    length(present) <= 3 || error("Relative permeability setup called for system with more than three phases.")
    has_ow = present.water && present.oil
    has_og = present.oil && present.gas
    has_wg = present.water && present.gas && !present.oil

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
        if has_ow
            krw[reg], krow[reg] = afi_relperm_pair(vals, :ow)
        end
        if has_og
            krg[reg], krog[reg] = afi_relperm_pair(vals, :og)
        end
        if has_wg
            krw[reg], krg[reg] = afi_relperm_pair(vals, :wg)
        end
    end
    w = remap_to_tuple(krw, kr_regs)
    g = remap_to_tuple(krg, kr_regs)
    ow = remap_to_tuple(krow, kr_regs)
    og = remap_to_tuple(krog, kr_regs)

    hyst = find_records(d, "KilloughRelPermHysteresis", "IX", steps = false, model = true, once = true)
    if !isnothing(hyst) && !disable_hysteresis
        h = hyst.value
        # Drainage, imbibition, hysteresis?
        tol = get(h, "ModificationParameter", 0.1)
        nw_hyst = JutulDarcy.KilloughHysteresis(tol = tol)
        wetting_hyst = String(get(h, "WettingPhaseHysteresis", "DRAINAGE"))
        if wetting_hyst == "DRAINAGE"
            w_hyst = JutulDarcy.NoHysteresis()
        elseif wetting_hyst == "IMBIBITION"
            w_hyst = JutulDarcy.ImbibitionOnlyHysteresis()
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
    if isnothing(rock_mgr) || disable_endscale
        hscale = false
        vscale = false
    else
        hscale = get(rock_mgr.value, "HorizontalEndPointScaling", false)
        vscale = get(rock_mgr.value, "VerticalEndPointScaling", false)
    end
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

function afi_relperm_pair(satfun, type::Symbol)
    if type == :ow
        tab_label_w = "WaterRelPermFunction"
        tab_label_nw = "OilInWaterRelPermFunction"
    elseif type == :og
        tab_label_w = "GasRelPermFunction"
        tab_label_nw = "OilInGasRelPermFunction"
    elseif type == :wg
        tab_label_w = "WaterRelPermFunction"
        tab_label_nw = "GasRelPermFunction"
    else
        error("Unsupported relperm pair type: $type")
    end
    if haskey(satfun, "RelPerm")
        kr_w = get_relperm_table(satfun, tab_label_w)
        kr_nw = get_relperm_table(satfun, tab_label_nw)
    elseif haskey(satfun, "CoreyRelPerm")
        tab_w = table_for_relperm(satfun, tab_label_w, "CoreyRelPerm")
        tab_nw = table_for_relperm(satfun, tab_label_nw, "CoreyRelPerm")
        kr_w, kr_nw = relperm_for_corey_pair(tab_w, tab_nw)
    end

    return (kr_w, kr_nw)
end

function get_relperm_table(satfun, name)
    krtab = table_for_relperm(satfun, name, "RelPerm")
    kr = krtab["RelPerm"]
    s = krtab["Saturation"]
    phase_label = krtab["Label"]
    return PhaseRelativePermeability(s, kr, label = phase_label)
end

function table_for_relperm(satfun, name, tabtype = "RelPerm")
    if name == "WaterRelPermFunction"
        phase_label = :w
    elseif name == "OilInWaterRelPermFunction"
        phase_label = :ow
    elseif name == "GasRelPermFunction"
        phase_label = :g
    elseif name == "OilInGasRelPermFunction"
        phase_label = :og
    else
        error("Unsupported relperm table name: $name")
    end
    krtab = get_saturation_function(satfun, tabtype, name, throw = true)
    krtab = copy(krtab)
    krtab["Label"] = phase_label
    return krtab
end

import JutulDarcy: brooks_corey_relperm

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
        if haskey(reservoir, Symbol(reg_no_label))
            reg_value = reservoir[Symbol(reg_no_label)]
        elseif length(keys(reg_map)) == 1
            jutul_message("Regions", "Only one region label found for region $reg_no_label, but did not find this region in the domain. Assuming all cells belong to this region.")
            nc = number_of_cells(reservoir)
            reg_value = ones(Int, nc)
        else
            error("Region $reg_no_label not found in reservoir domain.")
        end
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

function setup_corey_kr_afi(sr_tot, sr, n = 2.0, krmax = 1.0, krmax_end = krmax; label, npts = 1000)
    sr_other = sr_tot - sr
    ϵ = 1e-8
    ϵ = 0.0
    s = collect(range(sr, 1.0 - sr_other, length = npts))
    bc(s_i) = JutulDarcy.brooks_corey_relperm(s_i, n, sr + ϵ, krmax, sr_tot + 2*ϵ)
    kr = bc.(s)
    if krmax_end != krmax && sr_tot - sr > 0.0
        @assert krmax_end >= krmax
        push!(s, 1.0)
        push!(kr, krmax_end)
    end
    kr = PhaseRelativePermeability(s, kr, label = label)
    return kr
end

function residual_saturation(tab)
    return max(get(tab, "ResidualSaturation", 0.0), get(tab, "ConnateSaturation", 0.0))
end

function max_saturation(tab)
    return get(tab, "MaximumSaturation", 1.0)
end

function total_residual(tab, sr_tot)
    sr_from_max = 0.0
    return max(sr_tot, sr_from_max)
end

function table_krmax(tab)
    krmax = get(tab, "EndPointRelPerm", 1.0)
    krmax_crit = get(tab, "RelPermAtAssociatedCriticalSaturation", krmax)
    return (krmax, krmax_crit)
end

function relperm_for_corey_pair(tab_w, tab_nw)
    sr_w = residual_saturation(tab_w)
    sr_nw = residual_saturation(tab_nw)
    max_wetting = max_saturation(tab_w)

    n_w = get(tab_w, "Exponent", 2.0)
    n_nw = get(tab_nw, "Exponent", 2.0)

    sr_nw = max(sr_nw, 1.0 - max_wetting)
    sr_tot = sr_w + sr_nw

    sr_tot_w = total_residual(tab_w, sr_tot)
    sr_tot_nw = total_residual(tab_nw, sr_tot)

    krmax_w, krmax_crit_w = table_krmax(tab_w)
    krmax_nw, krmax_crit_nw = table_krmax(tab_nw)

    label_w = tab_w["Label"]
    label_nw = tab_nw["Label"]

    krw = setup_corey_kr_afi(sr_tot_w, sr_w, n_w, krmax_crit_w, krmax_w, label = label_w)
    krnw = setup_corey_kr_afi(sr_tot_nw, sr_nw, n_nw, krmax_crit_nw, krmax_nw, label = label_nw)
    return (krw, krnw)
end
