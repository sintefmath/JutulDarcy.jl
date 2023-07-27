function parse_keyword!(data, outer_data, units, f, ::Val{:WELSPECS})
    d = "Default"
    defaults = [d,     d,  -1,  -1, NaN,   d,     0.0, "STD", "SHUT", "YES",   0, "SEG", 0,     d, d, "STD"]
    utypes =   [:id, :id, :id, :id, :length, :id, :length,   :id,    :id,   :id, :id,   :id, :id, :id, :id, :id]
    data["WELSPECS"] = parse_defaulted_group(f, defaults)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:COMPDAT})
    d = "Default"
    defaults = [d, -1, -1, -1, -1, "OPEN", -1, -1.0, -1.0, -1.0, 0.0, -1, "Z", -1]
    compdat = parse_defaulted_group(f, defaults)
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

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONPROD})
    d = "Default"
    defaults = [d, "OPEN", d, Inf, Inf, Inf, Inf, NaN, 0.0, 0, 0.0]
    data["WCONPROD"] = parse_defaulted_group(f, defaults)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONHIST})
    d = "Default"
    defaults = [d, "OPEN", d, 0.0, 0.0, 0.0, 0, d, 0.0, 0.0, 0.0]
    data["WCONHIST"] = parse_defaulted_group(f, defaults)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONINJ})
    d = "Default"
    defaults = [d, d, "OPEN", d, Inf, Inf, 0.0, "NONE", -1, Inf, 0.0, 0.0]
    data["WCONINJ"] = parse_defaulted_group(f, defaults)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONINJE})
    d = "Default"
    defaults = [d, d, "OPEN", d, Inf, Inf, NaN, Inf, 0, 0.0, 0.0, 0, 0, 0]
    data["WCONINJE"] = parse_defaulted_group(f, defaults)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:TSTEP})
    data["TSTEP"] = parse_deck_vector(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:TUNING})
    skip_records(f, 3)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RPTSCHED})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:NUPCOL})
    rec = read_record(f)
    tdims = [3];
    data["NUPCOL"] = only(parse_defaulted_line(rec, tdims))
end

