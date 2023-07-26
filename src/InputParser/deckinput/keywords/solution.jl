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
    data["PRESSURE"] = parse_grid_vector(f, outer_data["GRID"]["cartDims"], Float64)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RS})
    data["RS"] = parse_grid_vector(f, outer_data["GRID"]["cartDims"], Float64)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:ACTNUM})
    data["ACTNUM"] = parse_grid_vector(f, get_cartdims(outer_data), Bool)
end

function parse_keyword!(data, outer_data, units, f, ::Val{:RSVD})
    data["RSVD"] = parse_deck_matrix(f)
end
