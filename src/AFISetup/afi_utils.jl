
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

function obsh_to_summary(obsh, t_seconds = missing; start_date = missing, smooth = false, alpha = 0.7, remove_bhp_missing = true)
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
    function set_well_outputs!(dest, shortname, longname, well)
        rat, tot = get_obsh_rate_and_cumulative(well, longname, obsh, t_s; smooth = smooth, alpha = alpha)
        if !ismissing(rat)
            dest["$(shortname)R"] = rat
            dest["$(shortname)T"] = tot
        end
        return dest
    end
    for well in keys(obsh["wells_interp"])
        swdata = Dict()
        set_well_outputs!(swdata, "WOP", "OIL_PRODUCTION", well)
        set_well_outputs!(swdata, "WWP", "WATER_PRODUCTION", well)
        set_well_outputs!(swdata, "WWI", "WATER_INJECTION", well)
        set_well_outputs!(swdata, "WLP", "LIQUID_PRODUCTION", well)
        set_well_outputs!(swdata, "WGP", "GAS_PRODUCTION", well)

        # Bottom hole pressure
        I_bhp = get(obsh["wells_interp"][well], "BOTTOM_HOLE_PRESSURE", missing)
        if !ismissing(I_bhp)
            if remove_bhp_missing
                I_bhp = remove_missing(I_bhp)
            end
            if smooth
                I_bhp = smooth_interpolant(I_bhp, alpha)
            end
            swdata["WBHP"] = I_bhp.(t_s)
        end
        summary_obs["VALUES"]["WELLS"][string(well)] = swdata
    end

    fld = Dict()
    function sum_wells!(dest, dest_name, k)
        val = missing
        for (i, w) in enumerate(keys(summary_obs["VALUES"]["WELLS"]))
            wval = get(summary_obs["VALUES"]["WELLS"][w], k, missing)
            if ismissing(wval)
                continue
            end
            wval = summary_obs["VALUES"]["WELLS"][w][k]
            if ismissing(val)
                val = copy(wval)
            else
                val .+= wval
            end
        end
        if !ismissing(val)
            dest[dest_name] = val
        end
        return dest
    end
    # Oil
    sum_wells!(fld, "FOPR", "WOPR")
    sum_wells!(fld, "FOPT", "WOPT")
    # Water
    sum_wells!(fld, "FWPR", "WWPR")
    sum_wells!(fld, "FWPT", "WWPT")
    sum_wells!(fld, "FWIR", "WWIR")
    sum_wells!(fld, "FWIT", "WWIT")
    # Gas
    sum_wells!(fld, "FGPR", "WGPR")
    sum_wells!(fld, "FGPT", "WGPT")
    # Liquid
    sum_wells!(fld, "FLPR", "WLPR")
    sum_wells!(fld, "FLPT", "WLPT")

    summary_obs["VALUES"]["FIELD"] = fld
    summary_obs["UNIT_SYSTEM"] = "si"
    return summary_obs
end

function get_obsh_rate_and_cumulative(well, prefix, obsh, t_s; smooth = false, alpha = 0.1)
    rate_key = "$(prefix)_RATE"
    cumulative_key = "$(prefix)_CUML"

    well_interp = obsh["wells_interp"][well]
    rate_interp = get(well_interp, rate_key, missing)
    if ismissing(rate_interp)
        rate = cumulative = missing
    else
        if smooth
            rate_interp = smooth_interpolant(rate_interp, alpha)
            # rate_val = rate_interp.F
            # rate_t = rate_interp.X
            # rate_vals_smooth = exponential_smoothing(rate_val, alpha)
            # rate_interp = get_linear_interpolant(rate_t, rate_vals_smooth)
        end
        rate = rate_interp.(t_s)

        recompute_cumulative = smooth || !haskey(well_interp, cumulative_key)
        if recompute_cumulative
            dt = diff([0.0; t_s])
            cumulative = cumsum(rate .* dt)
        else
            cumulative = well_interp[cumulative_key].(t_s)
        end
    end
    return (rate, cumulative)
end

function remove_missing(I::Jutul.LinearInterpolant)
    vals = I.F
    t = I.X
    new_vals = Float64[]
    new_t = Float64[]
    for (i, v) in enumerate(vals)
        if isfinite(v) && !isapprox(v, 0.0, atol = 1e-8)
            push!(new_vals, v)
            push!(new_t, t[i])
        end
    end
    if length(new_vals) <= 1
        return I
    else
        return Jutul.get_1d_interpolator(new_t, new_vals, cap_endpoints = true)
    end
end

function smooth_interpolant(I::Jutul.LinearInterpolant, alpha::Float64)
    vals = I.F
    smoothed_vals = exponential_smoothing(vals, alpha)
    return Jutul.get_1d_interpolator(I.X, smoothed_vals, cap_endpoints = true)
end

function exponential_smoothing(data::Vector{Float64}, alpha::Float64)
    alpha < 0.0 && error("Alpha must be in [0, 1], got $alpha")
    alpha > 1.0 && error("Alpha must be in [0, 1], got $alpha")
    smoothed = copy(data)
    for i in 2:length(data)
        smoothed[i] = alpha * data[i] + (1 - alpha) * smoothed[i - 1]
    end
    return smoothed
end
