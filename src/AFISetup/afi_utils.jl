
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

function obsh_to_summary(obsh, t_seconds = missing; start_date = missing)
    if ismissing(t_seconds)
        WI = obsh["wells_interp"]
        w_sample = first(values(WI))
        t_seconds = w_sample["seconds"]
    end
    summary_obs = Dict()
    summary_obs["TIME"] = (start_date = start_date, seconds = copy(t_seconds))
    summary_obs["VALUES"] = Dict()
    summary_obs["VALUES"]["WELLS"] = Dict()

    function get_data_obsh(name, t, cumulative = false)
        wdata = obsh["wells_interp"][name]
        t_i = t_jutul.*si_unit(:day)
        if t == :wwir
            if cumulative
                ix = to_cumulative(t_i, wdata["WATER_INJECTION_RATE"].(t_i))
            else
                ix = wdata["WATER_INJECTION_RATE"].(t_i).*si_unit(:day)
            end
        elseif t == :wwpr
            if cumulative
                ix = to_cumulative(t_i, wdata["WATER_PRODUCTION_RATE"].(t_i))
            else
                ix = wdata["WATER_PRODUCTION_RATE"].(t_i).*si_unit(:day)
            end
        elseif t == :wopr
            if cumulative
                ix = to_cumulative(t_i, wdata["OIL_PRODUCTION_RATE"].(t_i))
            else
                ix = wdata["OIL_PRODUCTION_RATE"].(t_i).*si_unit(:day)
            end
        elseif t == :lrat
            if cumulative
                ix = to_cumulative(t_i, wdata["LIQUID_PRODUCTION_RATE"].(t_i))
            else
                ix = wdata["LIQUID_PRODUCTION_RATE"].(t_i).*si_unit(:day)
            end
        elseif t == :wgpr
            return 0 .* t_i
            if cumulative
                ix = to_cumulative(t_i, wdata["GAS_PRODUCTION_RATE"].(t_i))
            else
                ix = wdata["GAS_PRODUCTION_RATE"].(t_i).*si_unit(:day)
            end
        elseif t == :wbhp
            if cumulative
                error("Cumulative BHP makes no sense")
            else
                ix = wdata["BOTTOM_HOLE_PRESSURE"].(t_i)./si_unit(:bar)
            end
        end
        return ix
    end

    t_s = summary_obs["TIME"].seconds
    dt_s = diff([0; t_s])
    for well in keys(obsh["wells_interp"])
        swdata = Dict()
        wopt = obsh["wells_interp"][well]["OIL_PRODUCTION_CUML"].(t_s)
        swdata["WOPT"] = wopt
        swdata["WOPR"] = diff([0; wopt])./dt_s
        wwpt = obsh["wells_interp"][well]["WATER_PRODUCTION_CUML"].(t_s)
        swdata["WWPT"] = wwpt
        swdata["WWPR"] = diff([0; wwpt])./dt_s
        wwit = obsh["wells_interp"][well]["WATER_INJECTION_CUML"].(t_s)
        swdata["WWIT"] = wwit
        swdata["WWIR"] = diff([0; wwit])./dt_s
        swdata["WBHP"] = obsh["wells_interp"][well]["BOTTOM_HOLE_PRESSURE"].(t_s)
        summary_obs["VALUES"]["WELLS"][string(well)] = swdata
    end
    summary_obs["UNIT_SYSTEM"] = "si"
    return summary_obs
end
