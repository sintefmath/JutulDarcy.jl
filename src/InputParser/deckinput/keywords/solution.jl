# Utilities

function get_cartdims(outer_data)
    g = get_section(outer_data, :GRID)
    @assert haskey(g, "cartDims") "Cannot access cartDims, has not been set."
    return g["cartDims"]
end

function set_cartdims!(outer_data, dim)
    @assert length(dim) == 3
    g = get_section(outer_data, :GRID)
    dim = tuple(dim...)
    gdata = get_section(outer_data, :GRID)
    gdata["cartDims"] = dim
    gdata["CURRENT_BOX"] = (lower = (1, 1, 1), upper = dim)
end

# Keywords follow

function parse_keyword!(data, outer_data, units, f, ::Val{:SGAS})
    data["SGAS"] = parse_grid_vector(f, outer_data["GRID"]["cartDims"], Float64)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:SWAT})
    data["SWAT"] = parse_grid_vector(f, outer_data["GRID"]["cartDims"], Float64)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:PRESSURE})
    p = parse_grid_vector(f, outer_data["GRID"]["cartDims"], Float64)
    swap_unit_system!(p, units, :pressure)
    data["PRESSURE"] = p
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RS})
    rs = parse_grid_vector(f, outer_data["GRID"]["cartDims"], Float64)
    swap_unit_system!(rs, units, :u_rs)
    data["RS"] = rs
end

function parse_keyword!(data, outer_data, units, f, ::Val{:ACTNUM})
    data["ACTNUM"] = parse_grid_vector(f, get_cartdims(outer_data), Bool)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RSVD})
    rs = parse_deck_matrix(f)
    swap_unit_system_axes!(rs, units, (:length, :u_rs))
    data["RSVD"] = rs
end

function parse_keyword!(data, outer_data, units, f, ::Val{:EQUIL})
    n = number_of_tables(outer_data, :equil)
    def = [0.0, NaN, 0.0, 0.0, 0.0, 0.0, 0, 0, 0]
    eunits = (:length, :pressure, :length, :pressure, :length, :pressure, :id, :id, :id)
    out = []
    for i = 1:n
        rec = read_record(f)
        result = parse_defaulted_line(rec, def)
        swap_unit_system_axes!(result, units, eunits)
        push!(out, result)
    end
    data["EQUIL"] = out
end
