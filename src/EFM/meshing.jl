fracture_mesh_options(fracture::AbstractEmbeddedFracture) = FractureMeshOptions()
fracture_mesh_options(fracture::Union{PolygonFracture, RectangleFracture, EllipticalFracture}) = fracture.mesh_options

function fracture_mesh_override(values, name, n)
    if values isa AbstractVector
        length(values) == n || throw(ArgumentError("$name must contain one entry per input fracture"))
        out = collect(values)
    else
        out = fill(values, n)
    end
    return out
end

function fracture_grids(fractures, frames, tol::T; resolution = nothing, cell_length = nothing) where T
    n = length(fractures)
    resolutions = fracture_mesh_override(resolution, "resolution", n)
    lengths = fracture_mesh_override(cell_length, "cell_length", n)
    grids = Vector{Union{Nothing, FractureGrid{T}}}(undef, n)
    for (id, fracture) in enumerate(fractures)
        r, h = resolutions[id], lengths[id]
        if isnothing(r) && isnothing(h)
            options = fracture_mesh_options(fracture)
            r, h = options.resolution, options.cell_length
        end
        options = FractureMeshOptions(; resolution = r, cell_length = h, float_type = T)
        grids[id] = fracture_grid(fracture, frames[id], options, tol)
    end
    return grids
end

function fracture_grid_axes(fracture::AbstractEmbeddedFracture, frame::FractureFrame{T}) where T
    return Point2{T}(1, 0), Point2{T}(0, 1)
end

function fracture_grid_axes(fracture::Union{RectangleFracture, EllipticalFracture}, frame::FractureFrame{T}) where T
    axis = fracture.axes[1]
    u = normalize(Point2{T}(dot(axis, frame.u), dot(axis, frame.v)))
    return u, Point2{T}(-u[2], u[1])
end

function fracture_grid(fracture, frame::FractureFrame{T}, options, tol) where T
    grid = nothing
    if !isnothing(options.resolution) || !isnothing(options.cell_length)
        u, v = fracture_grid_axes(fracture, frame)
        points = [Point2{T}(dot(p, u), dot(p, v)) for part in frame.parts for p in part]
        low = reduce((a, b) -> min.(a, b), points)
        high = reduce((a, b) -> max.(a, b), points)
        extent = high-low
        counts = options.resolution
        if isnothing(counts)
            counts = (max(1, ceil(Int, extent[1]/options.cell_length[1])),
                max(1, ceil(Int, extent[2]/options.cell_length[2])))
        end
        spacing = (extent[1]/counts[1], extent[2]/counts[2])
        minimum(spacing) > tol || throw(ArgumentError("fracture grid spacing must exceed the geometric tolerance; reduce resolution or use a smaller tolerance"))
        grid = FractureGrid(low[1]*u+low[2]*v, (u, v), spacing, counts)
    end
    return grid
end

function fracture_grid_coordinates(grid, p)
    delta = p-grid.origin
    return (dot(delta, grid.axes[1])/grid.spacing[1], dot(delta, grid.axes[2])/grid.spacing[2])
end

function fracture_grid_cell(grid, i, j)
    u, v = grid.axes
    a, b = grid.spacing[1]*u, grid.spacing[2]*v
    origin = grid.origin+(i-1)*a+(j-1)*b
    return [origin, origin+a, origin+a+b, origin+b]
end

function split_fracture_patch(patch, grid, frame, tol)
    T = typeof(tol)
    coordinates = [fracture_grid_coordinates(grid, p) for part in patch.parts for p in part]
    low = ntuple(d -> clamp(floor(Int, minimum(p -> p[d], coordinates))+1, 1, grid.resolution[d]), 2)
    high = ntuple(d -> clamp(floor(Int, maximum(p -> p[d], coordinates))+1, 1, grid.resolution[d]), 2)
    pieces = FracturePatch{T}[]
    tiles = Tuple{Int, Int}[]
    for j in low[2]:high[2], i in low[1]:high[1]
        boundary = fracture_grid_cell(grid, i, j)
        parts = Vector{Point2{T}}[]
        for part in patch.parts
            clipped = clip_polygon(part, boundary, tol)
            if length(clipped) >= 3 && polygon_area_centroid(clipped)[1] > tol^2
                push!(parts, clipped)
            end
        end
        for component in connected_patch_parts(parts, tol)
            area, center = patch_geometry(component, frame)
            origin = first(first(component))
            diameter = maximum(norm(p-origin) for part in component for p in part)
            if area > tol*diameter
                push!(pieces, FracturePatch(patch.fracture, copy(patch.matrix_cells), component,
                    center, area, patch.matrix_areas.*(area/patch.area), copy(patch.matrix_weights)))
                push!(tiles, (i, j))
            end
        end
    end
    return pieces, tiles
end

function fracture_component_root(parents, i)
    while parents[i] != i
        i = parents[i]
    end
    return i
end

function merged_fracture_patch(pieces, ids, frame, tol)
    T = typeof(tol)
    parts = Vector{Point2{T}}[]
    areas, storage = Dict{Int, T}(), Dict{Int, T}()
    for id in ids
        piece = pieces[id]
        append!(parts, piece.parts)
        for (k, cell) in enumerate(piece.matrix_cells)
            areas[cell] = get(areas, cell, zero(T))+piece.matrix_areas[k]
            storage[cell] = get(storage, cell, zero(T))+piece.area*piece.matrix_weights[k]
        end
    end
    area, center = patch_geometry(parts, frame)
    cells = sort!(collect(keys(areas)))
    weights = [storage[cell]/area for cell in cells]
    return FracturePatch(pieces[first(ids)].fracture, cells, parts, center, area,
        [areas[cell] for cell in cells], weights)
end

function fracture_grid_edge_ids(grid, i, j)
    nx, ny = grid.resolution
    edges = zeros(Int, 4)
    # Edge order matches fracture_grid_cell: bottom, right, top, left.
    if j > 1
        edges[1] = (nx-1)*ny+(i-1)*(ny-1)+j-1
    end
    if i < nx
        edges[2] = (i-1)*ny+j
    end
    if j < ny
        edges[3] = (nx-1)*ny+(i-1)*(ny-1)+j
    end
    if i > 1
        edges[4] = (i-2)*ny+j
    end
    return edges
end

function add_fracture_grid_edges!(face_patches, patch, id, tile, grid, frame::FractureFrame{T}, tol) where T
    i, j = tile
    vertices = fracture_grid_cell(grid, i, j)
    edges = fracture_grid_edge_ids(grid, i, j)
    for k in eachindex(edges)
        if edges[k] > 0
            a, b = vertices[k], vertices[mod1(k+1, 4)]
            segment_length = norm(b-a)
            direction = (b-a)/segment_length
            for (lo, hi) in line_intervals(patch.parts, a, direction, tol; limits = (0, segment_length))
                connected = get!(face_patches, (patch.fracture, -edges[k]),
                    Dict{Int, Vector{Tuple{Point3{T}, Point3{T}}}}())
                segments = get!(connected, id, Tuple{Point3{T}, Point3{T}}[])
                push!(segments, (unproject(frame, a+lo*direction), unproject(frame, a+hi*direction)))
            end
        end
    end
    return face_patches
end

function remesh_fracture_patches(patches::Vector{FracturePatch{T}}, face_patches, grids, frames, mesh, tol) where T
    if any(grid -> !isnothing(grid), grids)
        pieces = FracturePatch{T}[]
        tiles = Tuple{Int, Int}[]
        descendants = Vector{Int}[]
        for patch in patches
            start = length(pieces)+1
            grid = grids[patch.fracture]
            if isnothing(grid)
                push!(pieces, patch)
                push!(tiles, (0, 0))
            else
                split, split_tiles = split_fracture_patch(patch, grid, frames[patch.fracture], tol)
                append!(pieces, split)
                append!(tiles, split_tiles)
            end
            push!(descendants, collect(start:length(pieces)))
        end
        piece_faces = empty(face_patches)
        for ((fi, face), connected) in pairs(face_patches), (old, segments) in pairs(connected)
            tagged = [(face, a, b) for (a, b) in segments]
            for id in descendants[old]
                out = get!(piece_faces, (fi, face), empty(connected))
                dest = get!(out, id, Tuple{Point3{T}, Point3{T}}[])
                if isnothing(grids[fi])
                    append!(dest, segments)
                else
                    for (_, a, b) in clipped_face_segments(frames[fi], tagged, tol; parts = pieces[id].parts)
                        push!(dest, (a, b))
                    end
                end
            end
        end
        interfaces = fracture_interfaces(piece_faces, tol)
        parents = collect(eachindex(pieces))
        for ((l, r), segments) in sort!(collect(interfaces); by = first)
            if tiles[l] == tiles[r] && !isnothing(grids[pieces[l].fracture])
                a, b = fracture_component_root(parents, l), fracture_component_root(parents, r)
                parents[max(a, b)] = min(a, b)
            end
        end
        groups = Dict{Int, Vector{Int}}()
        for id in eachindex(pieces)
            push!(get!(groups, fracture_component_root(parents, id), Int[]), id)
        end
        patches = FracturePatch{T}[]
        renumber = zeros(Int, length(pieces))
        group_tiles = Tuple{Int, Int}[]
        for root in sort!(collect(keys(groups)))
            ids = groups[root]
            fi = pieces[root].fracture
            patch = pieces[root]
            if length(ids) > 1
                patch = merged_fracture_patch(pieces, ids, frames[fi], tol)
            end
            push!(patches, patch)
            renumber[ids] .= length(patches)
            push!(group_tiles, tiles[root])
        end
        face_patches = empty(piece_faces)
        for (key, connected) in pairs(piece_faces), (id, segments) in pairs(connected)
            out = get!(face_patches, key, empty(connected))
            append!(get!(out, renumber[id], Tuple{Point3{T}, Point3{T}}[]), segments)
        end
        for (id, patch) in enumerate(patches)
            grid = grids[patch.fracture]
            if !isnothing(grid)
                add_fracture_grid_edges!(face_patches, patch, id, group_tiles[id], grid, frames[patch.fracture], tol)
            end
        end
    end
    return patches, face_patches
end
