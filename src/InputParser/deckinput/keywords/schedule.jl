function parse_keyword!(data, outer_data, units, f, ::Val{:WELSPECS})
    d = "Default"
    defaults = [d,     d,  -1,  -1, NaN,   d,     0.0, "STD", "SHUT", "YES",   0, "SEG", 0,     d, d, "STD"]
    utypes =   [:id, :id, :id, :id, :length, :id, :length,   :id,    :id,   :id, :id,   :id, :id, :id, :id, :id]
    wspecs = parse_defaulted_group(f, defaults)
    for ws in wspecs
        swap_unit_system_axes!(ws, units, utypes)
    end
    data["WELSPECS"] = wspecs
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
    wconprod = parse_defaulted_group(f, defaults)
    utypes = identity_unit_vector(defaults)
    utypes[4] = :liquid_rate_surface
    utypes[5] = :liquid_rate_surface
    utypes[6] = :gas_rate_surface
    utypes[7] = :liquid_rate_surface
    utypes[8] = :liquid_rate_reservoir
    utypes[9] = :pressure
    utypes[10] = :pressure
    for wp in wconprod
        swap_unit_system_axes!(wp, units, utypes)
    end
    data["WCONPROD"] = wconprod
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONHIST})
    d = "Default"
    defaults = [d, "OPEN", d, 0.0, 0.0, 0.0, 0, d, 0.0, 0.0, 0.0]
    out = parse_defaulted_group(f, defaults)
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

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONINJ})
    d = "Default"
    defaults = [d, d, "OPEN", d, Inf, Inf, 0.0, "NONE", -1, Inf, 0.0, 0.0]
    out = parse_defaulted_group(f, defaults)
    utypes = identity_unit_vector(defaults)
    utypes[6] = :liquid_rate_surface
    utypes[9] = :pressure
    utypes[10] = :pressure
    for o in out
        # TODO: This should probably use the WELSPECS declaration
        if o[2] == "GAS"
            utypes[5] = :gas_rate_surface
            utypes[12] = :u_rv
        else
            utypes[5] = :liquid_rate_surface
            if o[2] == "OIL"
                utypes[12] = :u_rs
            end
        end
        swap_unit_system_axes!(o, units, utypes)
    end
    data["WCONINJ"] = out
end

function parse_keyword!(data, outer_data, units, f, ::Val{:WCONINJE})
    d = "Default"
    defaults = [d, d, "OPEN", d, Inf, Inf, NaN, Inf, 0, 0.0, 0.0, 0, 0, 0]
    out = parse_defaulted_group(f, defaults)
    utypes = identity_unit_vector(defaults)
    utypes[6] = :liquid_rate_reservoir
    utypes[7] = :pressure
    utypes[8] = :pressure
    for o in out
        # TODO: This should probably use the WELSPECS declaration
        if o[2] == "GAS"
            utypes[5] = :gas_rate_surface
            utypes[10] = :u_rv
        else
            utypes[5] = :liquid_rate_surface
            if o[2] == "OIL"
                utypes[12] = :u_rs
            end
        end
        swap_unit_system_axes!(o, units, utypes)
    end
    data["WCONINJE"] = out
end

function parse_keyword!(data, outer_data, units, f, ::Val{:TSTEP})
    dt = parse_deck_vector(f)
    swap_unit_system!(dt, units, :time)
    data["TSTEP"] = dt
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

