function finish_current_section!(data, units, cfg, outer_data, ::Val{:RUNSPEC})

end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NOECHO})
    # Do nothing
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ECHO})
    # Do nothing
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MESSAGES})
    parser_message(cfg, outer_data, "MESSAGES", PARSER_MISSING_SUPPORT)
    rec = read_record(f)
    # TODO: Process the record.
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:START})
    rec = read_record(f)
    tdims = [1, "JAN", 1970, "00:00:00"];
    start = parse_defaulted_line(rec, tdims)
    data["START"] = convert_date_kw(start)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TITLE})
    m = next_keyword!(f)
    data["TITLE"] = m
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:METRIC})
    data["METRIC"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FIELD})
    data["FIELD"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WATER})
    data["WATER"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:OIL})
    data["OIL"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:GAS})
    data["GAS"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:DISGAS})
    data["DISGAS"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:STONE1}, Val{:STONE2}, Val{:BAKER1}, Val{:BAKER2}})
    k = unpack_val(v)
    parser_message(cfg, outer_data, "$k", PARSER_JUTULDARCY_MISSING_SUPPORT)
    data["$k"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Val{:MISCSTR})
    parser_message(cfg, outer_data, "MISCSTR", PARSER_JUTULDARCY_MISSING_SUPPORT)
    rec = read_record(f)
    tdims = [NaN, NaN, NaN];
    l = parse_defaulted_line(rec, tdims)
    @assert !isnan(l[1])
    for i in 2:3
        if isnan(l[i])
            l[i] = l[1]
        end
    end
    cm = si_unit(:centi)*si_unit(:meter)
    d = si_unit(:dyne)
    for i in eachindex(l)
        l[i] *= d/cm
    end
    data["MISCSTR"] = l
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:AIM})
    data["AIM"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:HWELLS})
    data["HWELLS"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PRCORR})
    data["PRCORR"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MONITOR})
    data["MONITOR"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:THERMAL})
    data["THERMAL"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MECH})
    data["MECH"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:VAPOIL})
    data["VAPOIL"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:UNIFOUT})
    data["UNIFOUT"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:UNIFIN})
    data["UNIFIN"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:CPR})
    data["CPR"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NUMRES})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MULTSAVE})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TABDIMS})
    rec = read_record(f)
    tdims = [1, 1, 20, 20, 1, 20, 20, 1,
             1, -1, 10,  1, -1,  0,  0, -1,
             10, 10, 10, -1,  5,  5,  5,  0, -1];
    # TODO: Special logic for -1 entries
    data["TABDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:UDADIMS})
    rec = read_record(f)
    tdims = [0, 0, 100]
    data["UDADIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:UDQDIMS})
    rec = read_record(f)
    tdims = [16, 16, 0, 0, 0, 0, 0, 0, 0, 0, "N"]
    data["UDQDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:EXTRAPMS})
    rec = read_record(f)
    tdims = [0]
    data["EXTRAPMS"] = only(parse_defaulted_line(rec, tdims))
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FAULTDIM})
    rec = read_record(f)
    tdims = [0];
    data["FAULTDIM"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:GRIDOPTS})
    rec = read_record(f)
    tdims = ["NO", 0, 0];
    data["GRIDOPTS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ENDSCALE})
    rec = read_record(f)
    parser_message(cfg, outer_data, "ENDSCALE", PARSER_JUTULDARCY_MISSING_SUPPORT)
    tdims = ["NODIR", "REVERS", 1, 20, 0];
    data["ENDSCALE"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:EQLDIMS})
    rec = read_record(f)
    tdims = [1, 100, 50, 1, 50];
    data["EQLDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MSGFILE})
    rec = read_record(f)
    tdims = [1];
    data["MSGFILE"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:REGDIMS})
    rec = read_record(f)
    tdims = [1, 1, 0, 0, 0, 1, 0, 0, 0];
    data["REGDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WELLDIMS})
    rec = read_record(f)
    tdims = [0, 0, 0, 0, 5, 10, 5, 4, 3, 0, 1, 1, 10, 201]
    data["WELLDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:WSEGDIMS})
    rec = read_record(f)
    tdims = [0, 1, 1, 0]
    data["WSEGDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:VFPPDIMS})
    rec = read_record(f)
    tdims = [0, 0, 0, 0, 0, 0]
    data["VFPPDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:VFPIDIMS})
    rec = read_record(f)
    tdims = [0, 0, 0]
    data["VFPIDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:AQUDIMS})
    rec = read_record(f)
    tdims = [1, 1, 1, 36, 1, 1, 0, 0]
    data["AQUDIMS"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MISCIBLE})
    rec = read_record(f)
    tdims = [1, 20, "NONE"]
    parser_message(cfg, outer_data, "MISCIBLE", PARSER_JUTULDARCY_MISSING_SUPPORT)
    data["MISCIBLE"] = parse_defaulted_line(rec, tdims)
end
