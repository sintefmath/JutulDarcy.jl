function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NOECHO})
    # Do nothing
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ECHO})
    # Do nothing
end


function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MESSAGES})
    read_record(f)
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

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NONNC})
    data["NONNC"] = true
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

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TABDIMS})
    rec = read_record(f)
    tdims = [1, 1, 20, 20, 1, 20, 20, 1,
             1, -1, 10,  1, -1,  0,  0, -1,
             10, 10, 10, -1,  5,  5,  5,  0, -1];
    # TODO: Special logic for -1 entries
    data["TABDIMS"] = parse_defaulted_line(rec, tdims)
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

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:EQLDIMS})
    rec = read_record(f)
    tdims = [1, 100, 50, 1, 50];
    data["EQLDIMS"] = parse_defaulted_line(rec, tdims)
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
