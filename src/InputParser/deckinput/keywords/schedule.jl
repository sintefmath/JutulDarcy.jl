function finish_current_section!(data, units, cfg, outer_data, ::Val{:SCHEDULE})
    compord = outer_data["SCHEDULE"]["COMPORD"]
    for k in keys(get_wells(outer_data))
        if !haskey(compord, k)
            compord[k] = "TRACK"
        end
    end
end

function welspecs_to_well(ws)
    name = ws[1]::AbstractString
    # TODO: Parse more fields.
    w = (
        head = (ws[3], ws[4]),
        ref_depth = ws[5],
        preferred_phase = ws[6],
        drainage_radius = ws[7]
    )
    return (name, w)
end

function get_wells(outer_data)
    return outer_data["SCHEDULE"]["WELSPECS"]
end

function get_well(outer_data, name)
    return get_wells(outer_data)[name]
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WELSPECS})
    d = "Default"
    defaults = [d,     d,  -1,  -1, NaN,       d,     0.0, "STD", "SHUT", "YES",   0, "SEG",    0,  d,  d, "STD", 0]
    utypes =   [:id, :id, :id, :id, :length, :id, :length,   :id,    :id,   :id, :id,   :id, :id, :id, :id, :id, :id]
    wspecs = parse_defaulted_group(f, defaults)
    wells = get_wells(outer_data)
    for ws in wspecs
        swap_unit_system_axes!(ws, units, utypes)
        name, well = welspecs_to_well(ws)
        # @assert !haskey(wells, name) "Well $name is already defined."
        wells[name] = well
    end
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COMPORD})
    d = "Default"
    defaults = [d, "TRACK"]
    compord = parse_defaulted_group(f, defaults)
    out = outer_data["SCHEDULE"]["COMPORD"]
    for co in compord
        well, val = co
        out[well] = val
    end
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COMPDAT})
    d = "Default"
    defaults = [d, -1, -1, -1, -1, "OPEN", -1, NaN, NaN, NaN, 0.0, -1.0, "Z", -1.0]
    wells = get_wells(outer_data)
    compdat = parse_defaulted_group_well(f, defaults, wells, 1)
    # Unit conversion
    utypes = identity_unit_vector(defaults)
    utypes[8] = :transmissibility
    utypes[9] = :length
    utypes[10] = :Kh
    utypes[12] = :time_over_volume
    utypes[14] = :length
    for cd in compdat
        swap_unit_system_axes!(cd, units, utypes)
    end
    data["COMPDAT"] = compdat
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WELOPEN})
    d = "Default"
    defaults = [d, "OPEN", -1, -1, -1, -1, -1]
    wells = get_wells(outer_data)
    welopen = parse_defaulted_group_well(f, defaults, wells, 1)
    parser_message(cfg, outer_data, "WELTARG", PARSER_JUTULDARCY_MISSING_SUPPORT)
    data["WELOPEN"] = welopen
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WCONPROD})
    d = "Default"
    defaults = [
        d, "OPEN", d, -1.0, -1.0,
        -1.0, -1.0, -1.0, -1.0, -1.0,
        0, 0.0, NaN, NaN, NaN,
        NaN, NaN, NaN, NaN, NaN
        ]
    wells = get_wells(outer_data)
    wconprod = parse_defaulted_group_well(f, defaults, wells, 1)
    for kw in wconprod
        for i in 4:10
            # Handle defaults
            if kw[i] < 0
                kw[i] = Inf
            end
        end
    end
    utypes = identity_unit_vector(defaults)
    utypes[4] = :liquid_rate_surface
    utypes[5] = :liquid_rate_surface
    utypes[6] = :gas_rate_surface
    utypes[7] = :liquid_rate_surface
    utypes[8] = :liquid_rate_reservoir
    utypes[9] = :pressure
    for wp in wconprod
        swap_unit_system_axes!(wp, units, utypes)
    end
    data["WCONPROD"] = wconprod
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WCONHIST})
    d = "Default"
    defaults = [d, "OPEN", d, 0.0, 0.0, 0.0, 0, d, 0.0, 0.0, 0.0]
    wells = get_wells(outer_data)
    out = parse_defaulted_group_well(f, defaults, wells, 1)
    utypes = identity_unit_vector(defaults)
    utypes[4] = :liquid_rate_surface
    utypes[5] = :liquid_rate_surface
    utypes[6] = :gas_rate_surface
    utypes[9] = :pressure
    utypes[10] = :pressure
    utypes[11] = :gas_rate_surface
    for o in out
        swap_unit_system_axes!(o, units, utypes)
    end
    data["WCONHIST"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WCONINJ})
    d = "Default"
    defaults = [d, d, "OPEN", d, Inf, Inf, 0.0, "NONE", -1, Inf, 0.0, 0.0]
    wells = get_wells(outer_data)
    out = parse_defaulted_group_well(f, defaults, wells, 1)
    utypes = identity_unit_vector(defaults)
    utypes[6] = :liquid_rate_surface
    utypes[9] = :pressure
    utypes[10] = :pressure
    for o in out
        w = get_well(outer_data, o[1])
        if w.preferred_phase == "GAS"
            utypes[5] = :gas_rate_surface
            utypes[12] = :u_rv
        else
            utypes[5] = :liquid_rate_surface
            if w.preferred_phase == "OIL"
                utypes[12] = :u_rs
            end
        end
        swap_unit_system_axes!(o, units, utypes)
    end
    data["WCONINJ"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WCONINJE})
    d = "Default"
    defaults = [d, d, "OPEN", d, -1.0, -1.0, NaN, -1.0, 0, 0.0, 0.0, 0, 0, 0]
    wells = get_wells(outer_data)
    out = parse_defaulted_group_well(f, defaults, wells, 1)
    utypes = identity_unit_vector(defaults)
    utypes[6] = :liquid_rate_reservoir
    utypes[7] = :pressure
    utypes[8] = :pressure
    for o in out
        for i in [5, 6, 8]
            if o[i] <= 0
                o[i] = Inf
            end
        end
        w = get_well(outer_data, o[1])
        if w.preferred_phase == "GAS"
            utypes[5] = :gas_rate_surface
            utypes[10] = :u_rv
        else
            utypes[5] = :liquid_rate_surface
            if w.preferred_phase == "OIL"
                utypes[12] = :u_rs
            end
        end
        swap_unit_system_axes!(o, units, utypes)
    end
    data["WCONINJE"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WCONINJH})
    d = "Default"
    defaults = [d, d, "OPEN", 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, "RATE"]
    wells = get_wells(outer_data)
    out = parse_defaulted_group_well(f, defaults, wells, 1)
    utypes = identity_unit_vector(defaults)
    utypes[5] = :pressure
    utypes[6] = :pressure
    for o in out
        w = get_well(outer_data, o[1])
        if w.preferred_phase == "GAS"
            utypes[4] = :gas_rate_surface
            utypes[8] = :u_rv
        else
            utypes[4] = :liquid_rate_surface
            if w.preferred_phase == "OIL"
                utypes[8] = :u_rs
            end
        end
        swap_unit_system_axes!(o, units, utypes)
    end
    data["WCONINJH"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Union{Val{:WELLTARG}, Val{:WELTARG}})
    defaults = ["Default", "ORAT", NaN]
    wells = get_wells(outer_data)
    out = parse_defaulted_group_well(f, defaults, wells, 1)
    parser_message(cfg, outer_data, "WELTARG", PARSER_JUTULDARCY_MISSING_SUPPORT)
    data["WELTARG"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WEFAC})
    parser_message(cfg, outer_data, "WEFAC", PARSER_JUTULDARCY_PARTIAL_SUPPORT)

    defaults = ["Default", 1.0, "YES"]
    wells = get_wells(outer_data)
    parsed = parse_defaulted_group_well(f, defaults, wells, 1)
    push_and_create!(data, "WEFAC", parsed)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WPIMULT})
    parser_message(cfg, outer_data, "WPIMULT", PARSER_JUTULDARCY_MISSING_SUPPORT)

    defaults = ["Default", 1.0, -1, -1, -1, -1, -1, -1]
    wells = get_wells(outer_data)
    parsed = parse_defaulted_group_well(f, defaults, wells, 1)
    push_and_create!(data, "WPIMULT", parsed)
end

function convert_date_kw(t)
    @assert length(t) == 4
    function get_month(s)
        s = uppercase(s)
        if s == "JAN"
            return 1
        elseif s == "FEB"
            return 2
        elseif s == "MAR"
            return 3
        elseif s == "APR"
            return 4
        elseif s == "MAY"
            return 5
        elseif s == "JUN"
            return 6
        elseif s == "JLY" || s == "JUL"
            return 7
        elseif s == "AUG"
            return 8
        elseif s == "SEP"
            return 9
        elseif s == "OCT"
            return 10
        elseif s == "NOV"
            return 11
        else
            @assert s == "DEC" "Did not understand month format $s"
            return 12
        end
    end
    yr = t[3]
    m = get_month(t[2])
    d = t[1]

    return parse(DateTime, "$yr-$d-$m $(t[4])", dateformat"y-d-m H:M:S")
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:DATES})
    defaults = [1, "JAN", 1970, "00:00:00"]
    dt = parse_defaulted_group(f, defaults)
    out = Vector{DateTime}()
    sizehint!(out, length(dt))
    for t in dt
        if length(t[end]) > 8
            # TODO: Fix dirty hack for skipping ms
            t[end] = t[end][1:8]
        end
        push!(out, convert_date_kw(t))
    end
    data["DATES"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TSTEP})
    dt = parse_deck_vector(f)
    swap_unit_system!(dt, units, :time)
    data["TSTEP"] = dt
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TUNING})
    skip_records(f, 3)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTSCHED})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NUPCOL})
    rec = read_record(f)
    tdims = [3];
    data["NUPCOL"] = only(parse_defaulted_line(rec, tdims))
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WELLSTRE})
    n = compositional_number_of_components(outer_data)
    defaults = ["DEFAULT", 0.0]
    for i in 1:(n-1)
        push!(defaults, 0.0)
    end
    rec = parse_defaulted_group(f, defaults)
    for r in rec
        sum = 0.0
        for i in 1:n
            sum += r[i+1]
        end
        @assert sum ≈ 1.0 "Compositions in well stream $(r[1]) defined by WELLSTRE must sum up to one (was $sum)"
    end
    data["WELLSTRE"] = rec
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WINJGAS})
    defaults = ["Default", "GRUP", "Default", "Default", 0]
    wells = get_wells(outer_data)
    d = parse_defaulted_group_well(f, defaults, wells, 1)
    data["WINJGAS"] = d
end


function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WRFTPLT})
    parser_message(cfg, outer_data, "WRFTPLT", PARSER_MISSING_SUPPORT)
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:GRUPTREE})
    parser_message(cfg, outer_data, "GRUPTREE", PARSER_MISSING_SUPPORT)
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:AQUFETP})
    parser_message(cfg, outer_data, "AQUFETP", PARSER_MISSING_SUPPORT)
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:AQUANCON})
    parser_message(cfg, outer_data, "AQUANCON", PARSER_MISSING_SUPPORT)
    skip_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WELSEGS})
    d = "Default"
    # TODO: Last two entries for heat capacity / thermal conductivity are not
    # properly handled w.r.t. units.
    defaults = [d, NaN, 0.0, 1e-5, d, "HFA", "HO", 0.0, 0.0, 0.0, 0.0, 0.0]
    utypes = [:id, :length, :length, :volume, :id, :id, :id, :length, :length, :area, :id, :id]
    wells = get_wells(outer_data)
    rec = read_record(f)
    wheader = parse_defaulted_line(rec, defaults)
    swap_unit_system_axes!(wheader, units, utypes)

    wname = wheader[1]
    @assert wname != d
    wrec = wheader[2:end]

    # Now follows the segments
    defaults = []
    utypes = Symbol[]

    function add_entry!(v, u = :id)
        push!(defaults, v)
        push!(utypes, u)
    end

    add_entry!(-1)
    add_entry!(-1)
    add_entry!(-1)
    add_entry!(-1)
    add_entry!(0.0, :length)
    add_entry!(0.0, :length)
    add_entry!(0.0, :length)
    add_entry!(0.0, :area)
    add_entry!(0.0, :volume)
    add_entry!(0.0, :length)
    add_entry!(0.0, :length)
    add_entry!(0.0, :length)
    add_entry!(0.0, :area)
    add_entry!(0.0, :id)
    add_entry!(0.0, :id)

    segments = []
    while true
        rec = read_record(f)
        if length(rec) == 0
            break
        end
        l = parse_defaulted_line(rec, defaults)
        swap_unit_system_axes!(l, units, utypes)
        push!(segments, l)
    end
    if !haskey(data, "WELSEGS")
        data["WELSEGS"] = Dict{String, Any}()
    end
    data["WELSEGS"][wname] = (header = wheader, segments = segments)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COMPSEGS})
    rec = read_record(f)
    wname = only(parse_defaulted_line(rec, ["Default"]))
    @assert wname != "Default"

    defaults = []
    utypes = Symbol[]

    function add_entry!(v, u = :id)
        push!(defaults, v)
        push!(utypes, u)
    end

    add_entry!(-1)
    add_entry!(-1)
    add_entry!(-1)
    add_entry!(-1)
    add_entry!(0.0, :length)
    add_entry!(0.0, :length)
    add_entry!("Default")
    add_entry!(-1)
    add_entry!(0.0, :length)
    add_entry!(0.0, :length)
    add_entry!(0)

    segments = []
    while true
        rec = read_record(f)
        if length(rec) == 0
            break
        end
        l = parse_defaulted_line(rec, defaults)
        swap_unit_system_axes!(l, units, utypes)
        push!(segments, l)
    end
    if !haskey(data, "COMPSEGS")
        data["COMPSEGS"] = Dict{String, Any}()
    end
    data["COMPSEGS"][wname] = segments
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:UDQ})
    skip_record(f)
    parser_message(cfg, outer_data, "UDQ", PARSER_MISSING_SUPPORT)
end
