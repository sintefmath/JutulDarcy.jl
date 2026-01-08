
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

    t_s = t_seconds
    for well in keys(obsh["wells_interp"])
        swdata = Dict()
        I = obsh["wells_interp"][well]
        swdata["WOPT"] = I["OIL_PRODUCTION_CUML"].(t_s)
        swdata["WOPR"] = I["OIL_PRODUCTION_RATE"].(t_s)

        swdata["WWPT"] = I["WATER_PRODUCTION_CUML"].(t_s)
        swdata["WWPR"] = I["WATER_PRODUCTION_RATE"].(t_s)

        swdata["WWIT"] = I["WATER_INJECTION_CUML"].(t_s)
        swdata["WWIR"] = I["WATER_INJECTION_RATE"].(t_s)
        swdata["WBHP"] = I["BOTTOM_HOLE_PRESSURE"].(t_s)

        swdata["WLPR"] = I["LIQUID_PRODUCTION_RATE"].(t_s)
        swdata["WLPT"] = I["LIQUID_PRODUCTION_CUML"].(t_s)

        summary_obs["VALUES"]["WELLS"][string(well)] = swdata
    end
    summary_obs["UNIT_SYSTEM"] = "si"
    return summary_obs
end
