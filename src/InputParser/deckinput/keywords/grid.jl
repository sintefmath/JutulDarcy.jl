
function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:INIT})
    data["INIT"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COORDSYS})
    read_record(f)
    @warn "COORDSYS skipped."
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COORD})
    coord = parse_deck_vector(f, Float64)
    coord = swap_unit_system!(coord, units, Val(:length))
    data["COORD"] = coord
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ZCORN})
    zcorn = parse_deck_vector(f, Float64)
    zcorn = swap_unit_system!(zcorn, units, Val(:length))
    data["ZCORN"] = zcorn
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:PERMX}, Val{:PERMY}, Val{:PERMZ}})
    k = unpack_val(v)
    vals = parse_grid_vector(f, get_cartdims(outer_data), Float64)
    vals = swap_unit_system!(vals, units, Val(:permeability))
    data["$k"] = vals
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PORO})
    data["PORO"] = parse_grid_vector(f, get_cartdims(outer_data), Float64)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NTG})
    data["NTG"] = parse_grid_vector(f, get_cartdims(outer_data), Float64)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:DX}, Val{:DY}, Val{:DZ}})
    k = unpack_val(v)
    data["$k"] = parse_grid_vector(f, get_cartdims(outer_data), Float64)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TOPS})
    tops = parse_deck_vector(f, Float64)
    data["TOPS"] = swap_unit_system!(tops, units, Val(:length))
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:DIMENS})
    rec = read_record(f)
    to_int = x -> Parsers.parse(Int, x)
    d = to_int.(filter!(x -> length(x) > 0, split(only(rec), DECK_SPLIT_REGEX)))
    data["DIMENS"] = d
    set_cartdims!(outer_data, d)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:SPECGRID})
    rec = read_record(f)
    tdims = [1, 1, 1, 1, "F"]
    data["SPECGRID"] = parse_defaulted_line(rec, tdims)
    set_cartdims!(outer_data, data["SPECGRID"][1:3])
end

