
function parse_keyword!(data, outer_data, units, f, ::Val{:RPTPROPS})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, f, v::Union{Val{:SWOF}, Val{:SGOF}})
    k = unpack_val(v)
    sat_tab = parse_saturation_table(f, outer_data)
    for tab in sat_tab
        swap_unit_system_axes!(tab, units, (:identity, :identity, :identity, :pressure))
    end
    data["$k"] = sat_tab
end


function parse_keyword!(data, outer_data, units, f, ::Val{:PVDG})
    pvdg = parse_dead_pvt_table(f, outer_data)
    for tab in pvdg
        swap_unit_system_axes!(tab, units, (:pressure, :gas_formation_volume_factor, :viscosity))
    end
    data["PVDG"] = pvdg
end

function parse_keyword!(data, outer_data, units, f, ::Val{:PVTO})
    data["PVTO"] = parse_live_pvt_table(f, outer_data)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:PVTW})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN]
    # TODO: This should handle regions
    data["PVTW"] = parse_defaulted_line(rec, tdims)
    @assert all(isfinite, data["PVTW"]) "PVTW cannot be defaulted"
end

function parse_keyword!(data, outer_data, units, f, ::Val{:PVCDO})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN]
    # TODO: This should handle regions, merge with PVTW etc
    data["PVCDO"] = parse_defaulted_line(rec, tdims)
    @assert all(isfinite, data["PVCDO"]) "PVCDO cannot be defaulted"
end

function parse_keyword!(data, outer_data, units, f, ::Val{:ROCK})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN, NaN]
    data["ROCK"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:DENSITY})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN]
    data["DENSITY"] = parse_defaulted_line(rec, tdims)
end
