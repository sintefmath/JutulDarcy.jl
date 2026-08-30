function table_to_relperm(swof; swcon = 0.0, first_label = :w, second_label = :ow)
    sw = vec(swof[:, 1])
    krw = vec(swof[:, 2])
    krw = PhaseRelativePermeability(sw, krw, label = first_label)
    kro = vec(swof[end:-1:1, 3])
    so = 1 .- sw
    so = vec(so[end:-1:1])
    @. so = so - swcon
    krow = PhaseRelativePermeability(so, kro, label = second_label)
    return (krw, krow)
end

function saturation_table_handle_defaults(s, f)
    if any(isnan, f)
        # NaN values are removed due to INPUT file shenanigans
        ix = findall(!isnan, f)
        s = s[ix]
        f = f[ix]
    end
    return (s, f)
end

function add_missing_endpoints(s, kr)
    copied = false
    if s[1] > 0.0
        copied = true
        s = vcat(0.0, s)
        kr = vcat(0.0, kr)
    end
    if s[end] < 1.0
        copied = true
        s = vcat(s, 1.0)
        kr = vcat(kr, kr[end])
    end
    if !copied
        s = copy(s)
        kr = copy(kr)
    end
    return (s, kr)
end

function ensure_endpoints!(x, f, ϵ)
    n = length(x)
    for i in (n-1):-1:2
        if f[i] != f[i-1]
            x[i] -= ϵ
            break
        end
    end
    for i in 1:(n-1)
        if f[i] != f[i+1]
            x[i] -= ϵ
            break
        end
    end
end

"""
    summary_result(case::JutulCase, res::ReservoirSimResult) # units from case input data
    summary_result(case::JutulCase, res::ReservoirSimResult, :field) # field units

Create a summary-like result as a Dict. This can subsequently be written to disk
using `GeoEnergyIO.write_jutuldarcy_summary`. The `usys` argument is used to
specify the unit system to use. If `usys` is `missing`, the unit system is
chosen based on the input data as either :field, :lab or :metric, if data is
present, otherwise it will be set to :metric.

# Examples
```julia
smry_jutul = summary_result(case, res, :field)
GeoEnergyIO.write_jutuldarcy_summary("FILENAME", smry_jutul, unified = true)
```
"""
function summary_result(case::Jutul.JutulCase, res::ReservoirSimResult, usys = missing; kwarg...)
    return summary_result(case.model, res.wells, res.states, usys; input_data = case.input_data, kwarg...)
end

"""
    summary_result(model::MultiModel, wellresult, states)

Low level version of `summary_result` that works directly on a `MultiModel`,
"""
function summary_result(model::MultiModel, wellresult, states = missing, usys = missing;
        input_data = missing,
        wells = !ismissing(wellresult),
        start_date = missing,
        field = wells
    )
    function to_summary(x::Dict)
        x_c = Dict{String, Vector{Float64}}()
        # TODO: Unit conversion here.
        for (k, v) in pairs(x)
            if k == :time
                continue
            end
            skey = uppercase(string(k))
            x_c[skey] = v.values
        end
        return x_c
    end

    data = input_data
    has_data = !isnothing(data) && !ismissing(data) && haskey(data, "RUNSPEC") 
    if ismissing(usys)
        if has_data
            rs = data["RUNSPEC"]
            if haskey(rs, "FIELD")
                usys = :field
            elseif haskey(rs, "LAB")
                usys = :lab
            elseif haskey(rs, "PVT-M")
                error("PVT-M not supported for unit conversion.")
            else
                usys = :metric
            end
        else
            usys = :metric
        end
    end

    out = Dict()
    out["UNIT_SYSTEM"] = string(usys)
    out["VALUES"] = vals = Dict()
    function get_values(t; kwarg...)
            rm = JutulDarcy.reservoir_measurables(
            model, wellresult, states;
            units = usys,
            type = t,
            kwarg...
        )
        return rm
    end
    t = missing
    if field
        f_smry = get_values(:field)
        vals["FIELD"] = to_summary(f_smry)
        t = f_smry[:time]
    end
    if wells
        w_smry = Dict()
        for w in keys(wellresult.wells)
            wi = get_values(:well, wells = w)
            w_smry["$w"] = to_summary(wi)
            if ismissing(t)
                t = wi[:time]
            end
        end
        vals["WELLS"] = w_smry
    end

    if ismissing(start_date)
        if has_data && haskey(data["RUNSPEC"], "START")
            start_date = data["RUNSPEC"]["START"]
        else
            start_date = missing
        end
    end
    out["TIME"] = (start_date = start_date, seconds = t)

    reservoir = reservoir_domain(model)
    mesh = physical_representation(reservoir)
    if mesh isa UnstructuredMesh && mesh.structure isa CartesianIndex
        dims = mesh.structure.I
    else
        dims = (number_of_cells(mesh), 1, 1)
    end
    out["DIMENS"] = dims
    return out
end

"""
    s = convert_summary(s_flat_fmt)
    s = convert_summary(s_flat_fmt, unit_system = "METRIC")

Convert a summary dictionary in flat format (as produced by
`GeoEnergy.read_summary`) to the JutulDarcy style summary format that includes
nested dictionaries for fields, and construct a `TIME` entry in seconds from the
`YEARS` vector. The original unit system indicated by `UNIT_SYSTEM` is preserved
for all value arrays. This is then suitable for writing to disk using
`GeoEnergyIO.write_jutuldarcy_summary` or plotting in `plot_summary`.
"""
function convert_summary(s; unit_system = get(s, "UNIT_SYSTEM", "METRIC"))
    if haskey(s, "TIME")
        # Already in Jutul style format.
        return s
    end
    out = Dict{String, Any}()
    out["VALUES"] = Dict{String, Any}(
        "FIELD" => Dict{String, Any}(),
        "GROUP" => Dict{String, Any}(),
        "WELLS" => Dict{String, Any}(),
    )
    out["UNIT_SYSTEM"] = unit_system
    for (k, v) in pairs(s)
        if v isa AbstractArray
            out["VALUES"]["FIELD"][k] = v
        elseif v isa Dict
            for (wk, wv) in pairs(v)
                if startswith(k, "W")
                    if !haskey(out["VALUES"]["WELLS"], wk)
                        out["VALUES"]["WELLS"][wk] = Dict{String, Any}()
                    end
                    out["VALUES"]["WELLS"][wk][k] = wv
                elseif wv isa AbstractArray
                    out["VALUES"]["FIELD"][k] = v
                end
            end
        end
    end
    seconds = s["YEARS"].*si_unit(:year)
    out["TIME"] = (start_date = get(s, "START_DATE", nothing), seconds = seconds)
    return out
end

function convert_summary_from_data_file(data::AbstractDict)

    rs = data["RUNSPEC"]
    sched = data["SCHEDULE"]
    start_date = get(rs, "START", nothing)

    wells_scattered = Dict{String, Any}()
    for k in keys(sched["WELSPECS"])
        wells_scattered[k] = Dict(
            :time => Float64[],
            :qws => Float64[],
            :qos => Float64[],
            :qgs => Float64[],
            :bhp => Float64[]
        )
    end
    timesteps = Float64[]
    current_time = 0.0

    function strip_defaulted(x)
        if isfinite(x)
            return x
        else
            return 0.0
        end
    end

    function add_well_entry(qos::Real, qws::Real, qgs::Real, bhp::Real, wellname::AbstractString, sgn::Real)
        dest = wells_scattered[wellname]
        push!(dest[:time], current_time)
        push!(dest[:qos], sgn*strip_defaulted(qos))
        push!(dest[:qws], sgn*strip_defaulted(qws))
        push!(dest[:qgs], sgn*strip_defaulted(qgs))
        push!(dest[:bhp], strip_defaulted(bhp))
    end

    function add_well_entry(keywords::Vector, idx_qws::Int, idx_qos::Int, idx_qgs::Int, idx_bhp::Int, sgn::Real)
        for kword in keywords
            wellname = kword[1]
            qos = kword[idx_qos]
            qws = kword[idx_qws]
            qgs = kword[idx_qgs]
            bhp = kword[idx_bhp]
            add_well_entry(qos, qws, qgs, bhp, wellname, sgn)
        end
    end
    for step in sched["STEPS"]
        for (key, kword) in pairs(step)
            if key == "DATES"
                for date in kword
                    cdate = start_date + Second(current_time)
                    dt = Float64(Second(date - cdate).value)
                    push!(timesteps, dt)
                    current_time += dt
                end
            elseif key == "TIME"
                for time in kword
                    dt = time - current_time
                    push!(timesteps, dt)
                    current_time = time
                end
            elseif key == "TSTEP"
                found_time = true
                for dt in kword
                    push!(timesteps, dt)
                    current_time = current_time + dt
                end
            elseif key == "WCONHIST"
                # WCONHIST
                # 4 - qos
                # 5 - qws
                # 6 - qgs
                # 10 - bhp
                add_well_entry(kword, 5, 4, 6, 10, -1)
            elseif key == "WCONPROD"
                # WCONPROD
                # 4 - qos
                # 5 - qws
                # 6 - qgs
                # 7 - lrat
                # 8 - resv
                # 9 - bhp
                add_well_entry(kword, 5, 4, 6, 9, -1)
            elseif key == "WCONINJE"
                # WCONINJE
                # 2 - type
                # 5 - rate
                # 7 - bhp
                for kw in kword
                    rate = kw[5]
                    phase = kw[2]
                    if phase == "WATER"
                        qws = rate
                        qos = 0.0
                        qgs = 0.0
                    elseif phase == "OIL"
                        qws = 0.0
                        qos = rate
                        qgs = 0.0
                    elseif phase == "GAS"
                        qws = 0.0
                        qos = 0.0
                        qgs = rate
                    else
                        error("Unknown phase $phase in WCONINJE")
                    end
                    add_well_entry(qos, qws, qgs, kw[7], kw[1], 1)
                end
            elseif key == "WCONINJH"
                # TODO: handle WCONINJH
                error("WCONINJH not yet implemented in convert_summary_from_data_file")
            end
        end
    end
    # Now unify and output as actual summary
    wells = Dict{String, Any}()
    out = Dict{String, Any}()
    out["VALUES"] = Dict(
        "FIELD" => Dict{String, Any}(),
        "GROUP" => Dict{String, Any}(),
        "WELLS" => wells,
    )
    seconds = cumsum([0; timesteps])
    all(isfinite, seconds) || error("Non-finite time values in summary conversion.")
    function sample(well, response, sgn)
        ws = wells_scattered[well]
        t = ws[:time]

        i_u = unique(i -> t[i], eachindex(t))
        t = t[i_u]
        r = sgn.*ws[response][i_u]
        r = max.(r, 0.0)
        I = get_1d_interpolator(t, r)
        if length(t) < 2
            val = zeros(length(seconds))
        else
            val = I.(seconds)
        end
        return val
    end
    for wellname in keys(wells_scattered)
        well = Dict{String, Any}()
        well["WBHP"] = sample(wellname, :bhp, 1.0)
        # Production
        well["WOPR"] = sample(wellname, :qos, -1.0)
        well["WGPR"] = sample(wellname, :qgs, -1.0)
        well["WWPR"] = sample(wellname, :qws, -1.0)
        # Injection
        well["WWIR"] = sample(wellname, :qws, 1.0)
        well["WOIR"] = sample(wellname, :qos, 1.0)
        well["WGIR"] = sample(wellname, :qgs, 1.0)

        wells[wellname] = well
    end

    out["TIME"] = (start_date = start_date, seconds = seconds)
    out["UNIT_SYSTEM"] = "SI"
    return out
end

function expand_summary(summary)
    # Add missing fields
end
