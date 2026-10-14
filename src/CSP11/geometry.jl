const spe11_geometry_cache = Ref{Any}(nothing)

find_closest_point(points, reference) = last(findmin(point -> sum(abs2, point - reference), points))

"""
    spe11c_elevation_offset(y)

Elevation offset (metres) used by the SPE11C reference-to-physical mapping.
"""
function spe11c_elevation_offset(y)
    f = (y - 2500.0)/2500.0
    return 150.0*(1.0 - f*f) + y/500.0
end

function parse_int_list(s)
    return parse.(Int, strip.(split(s, ',')))
end

function read_spe11_geometry()
    pth = joinpath(artifact"spe11a_geo", "spe11a.geo")
    # pth = joinpath(@__DIR__, "..", "data", "spe11a.geo")
    points = Dict{Int, SVector{2, Float64}}()
    lines = Dict{Int, Tuple{Int, Int}}()
    loops = Dict{Int, Vector{Int}}()
    surface_facies = Dict{Int, Int}()

    for line in eachline(pth)
        if (m = match(r"^Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),", line)) !== nothing
            i = parse(Int, m.captures[1])
            points[i] = SVector(parse(Float64, m.captures[2]), parse(Float64, m.captures[3]))
        elseif (m = match(r"^Line\((\d+)\)\s*=\s*\{([^}]+)\}", line)) !== nothing
            i = parse(Int, m.captures[1])
            p = parse_int_list(m.captures[2])
            lines[i] = (p[1], p[2])
        elseif (m = match(r"^Curve Loop\((\d+)\)\s*=\s*\{([^}]+)\}", line)) !== nothing
            loops[parse(Int, m.captures[1])] = parse_int_list(m.captures[2])
        elseif (m = match(r"^Physical Surface\([^,]+,\s*(\d+)\)\s*=\s*\{([^}]+)\}", strip(line))) !== nothing
            facies = parse(Int, m.captures[1])
            for surface in parse_int_list(m.captures[2])
                surface_facies[surface] = facies
            end
        end
    end

    polygons = NamedTuple[]
    for (surface, edges) in sort!(collect(loops), by = first)
        ids = Int[]
        for signed_edge in edges
            a, b = lines[abs(signed_edge)]
            if signed_edge < 0
                a, b = b, a
            end
            if isempty(ids)
                append!(ids, (a, b))
            else
                ids[end] == a || error("Invalid SPE11 geometry: curve loop $surface is not continuous")
                push!(ids, b)
            end
        end
        first(ids) == last(ids) || error("Invalid SPE11 geometry: curve loop $surface is not closed")
        polygon = [points[i] for i in ids[1:end-1]]
        xs = first.(polygon)
        ys = last.(polygon)
        push!(polygons, (
            facies = surface_facies[surface],
            polygon = polygon,
            xmin = minimum(xs), xmax = maximum(xs),
            ymin = minimum(ys), ymax = maximum(ys)
        ))
    end
    return polygons
end

function spe11_geometry()
    if isnothing(spe11_geometry_cache[])
        spe11_geometry_cache[] = read_spe11_geometry()
    end
    return spe11_geometry_cache[]
end

"""
    spe11_reference_volumes(case)

Physical and effective storage volumes integrated from the bundled facies
polygons, including the lateral boundary buffers. Pore volumes include the
small facies-7 porosity used to keep the flow equations well-defined.
"""
function spe11_reference_volumes(case)
    check_spe11_case(case)
    areas = zeros(7)
    side_lengths = zeros(7)
    spec = spe11_case_specification(case)
    scale_x, scale_z = spec.geometry_scale
    xmax = 2.8
    for p in spe11_geometry()
        polygon = p.polygon
        signed_area = 0.0
        for i in eachindex(polygon)
            a, b = polygon[i], polygon[mod1(i + 1, length(polygon))]
            signed_area += a[1]*b[2] - b[1]*a[2]
            if (isapprox(a[1], 0.0; atol = 1e-10) && isapprox(b[1], 0.0; atol = 1e-10)) ||
                    (isapprox(a[1], xmax; atol = 1e-10) && isapprox(b[1], xmax; atol = 1e-10))
                side_lengths[p.facies] += scale_z*abs(b[2] - a[2])
            end
        end
        areas[p.facies] += abs(signed_area)*scale_x*scale_z/2
    end
    thickness = spec.physical_size[2]
    boundary_areas = side_lengths.*thickness
    if case == :c
        # The C mapping preserves both volume and the areas of the lateral
        # boundary planes. Front and back each contain the B cross-section.
        boundary_areas .+= 2 .* areas
    end
    facies = spec.buffer_facies
    _, porosity = rock_props_from_satnum(collect(1:7), case)
    physical_bulk = sum(areas)*thickness
    physical_pore = sum(areas.*porosity)*thickness
    if case == :a
        buffer_bulk = 0.0
        buffer_pore = 0.0
    else
        buffer_bulk = spe11_buffer_length*sum(boundary_areas[facies])
        buffer_pore = spe11_buffer_length*sum((boundary_areas.*porosity)[facies])
    end
    return (physical_bulk = physical_bulk, physical_pore = physical_pore,
        buffer_bulk = buffer_bulk, buffer_pore = buffer_pore,
        bulk = physical_bulk + buffer_bulk, pore = physical_pore + buffer_pore)
end

function point_on_segment(x, y, a, b; atol = 100*eps(Float64))
    cross = (x - a[1])*(b[2] - a[2]) - (y - a[2])*(b[1] - a[1])
    abs(cross) <= atol || return false
    return min(a[1], b[1]) - atol <= x <= max(a[1], b[1]) + atol &&
           min(a[2], b[2]) - atol <= y <= max(a[2], b[2]) + atol
end

function point_in_polygon(x, y, polygon)
    inside = false
    j = lastindex(polygon)
    for i in eachindex(polygon)
        pi = polygon[i]
        pj = polygon[j]
        point_on_segment(x, y, pi, pj) && return true
        if (pi[2] > y) != (pj[2] > y)
            x_intersection = (pj[1] - pi[1])*(y - pi[2])/(pj[2] - pi[2]) + pi[1]
            if x < x_intersection
                inside = !inside
            end
        end
        j = i
    end
    return inside
end

"""
    spe11_facies(x, elevation; case = :b)

Return the SPE11 facies number (1--7) at a point in the A or B reference geometry.
Coordinates are in metres. SPE11C points must first be mapped back to the
reference configuration (see the three-argument method).
"""
function spe11_facies(x, elevation; case = :b)
    check_spe11_case(case)
    scale_x, scale_z = spe11_case_specification(case).geometry_scale
    x_geo = x/scale_x
    z_geo = elevation/scale_z
    for p in spe11_geometry()
        if p.xmin <= x_geo <= p.xmax && p.ymin <= z_geo <= p.ymax &&
                point_in_polygon(x_geo, z_geo, p.polygon)
            return p.facies
        end
    end
    throw(ArgumentError("Point ($x, $elevation) lies outside the SPE11 geometry"))
end

"""
    spe11_facies(x, y, depth; case = :c)

Return the facies at a point expressed in Jutul's depth-positive coordinate
system. For case C, `depth` is physical depth and the arch mapping is inverted
before querying the two-dimensional reference geometry.
"""
function spe11_facies(x, y, depth; case = :c)
    check_spe11_case(case)
    if case == :c
        offset = spe11c_elevation_offset(y)
    else
        offset = 0.0
    end
    height = spe11_case_specification(case).physical_size[3]
    reference_elevation = height - depth - offset
    return spe11_facies(x, reference_elevation; case = case)
end
