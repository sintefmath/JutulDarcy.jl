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
    summary_result(case::JutulCase, res::ReservoirSimResult, usys = missing)

Write a summary-like result to a Dict. This can subsequently be written to disk
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
function summary_result(case::JutulCase, res::ReservoirSimResult, usys = missing)
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

    data = case.input_data
    has_data = !isnothing(data) && haskey(data, "RUNSPEC") 
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
    out["VALUES"] = vals = Dict()
    function get_values(t; kwarg...)
            rm = JutulDarcy.reservoir_measurables(
            case, res;
            units = usys,
            type = t,
            kwarg...
        )
        return rm
    end
    f_smry = get_values(:field)
    vals["FIELD"] = to_summary(f_smry)
    w_smry = Dict()
    for w in keys(res.wells.wells)
        wi = get_values(:well, wells = w)
        w_smry["$w"] = to_summary(wi)
    end
    vals["WELLS"] = w_smry

    if has_data && haskey(data["RUNSPEC"], "START")
        start = data["RUNSPEC"]["START"]
    else
        start = missing
    end
    out["TIME"] = (start_date = start, seconds = f_smry[:time])

    reservoir = reservoir_domain(case)
    mesh = physical_representation(reservoir)
    if mesh.structure isa CartesianIndex
        dims = mesh.structure.I
    else
        dims = (number_of_cells(mesh), 1, 1)
    end
    out["DIMENS"] = dims
    return out
end
