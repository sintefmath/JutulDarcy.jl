function fracture_frame(f::AbstractEmbeddedFracture, tol)
    T = typeof(tol)
    pts = Point3{T}.(fracture_polygon(f))
    origin, normal = first(pts), Point3{T}(f.normal)
    u, v = plane_basis(normal)
    parts = convex_parts([Point2{T}(dot(p-origin, u), dot(p-origin, v)) for p in pts], tol)
    return FractureFrame(origin, normal, u, v, parts)
end

function patch_geometry(parts, frame)
    T = eltype(frame.origin)
    area = zero(T)
    centroid = zero(Point2{T})
    for p in parts
        A, c = polygon_area_centroid(p)
        area += A
        centroid += A*c
    end
    return area, unproject(frame, centroid/area)
end

function clipped_face_segments(frame, segments, tol; parts = frame.parts)
    T = eltype(frame.origin)
    out = Tuple{Int, Point3{T}, Point3{T}}[]
    for (face, a, b) in segments
        pa, pb = project(frame, a), project(frame, b)
        L = norm(pb-pa)
        L > tol || continue
        direction = (pb-pa)/L
        for (lo, hi) in line_intervals(parts, pa, direction, tol; limits = (0, L))
            push!(out, (face, unproject(frame, pa+lo*direction), unproject(frame, pa+hi*direction)))
        end
    end
    return out
end

function connected_patch_parts(parts, tol)
    T = typeof(tol)
    length(parts) == 1 && return [parts]
    components = collect(eachindex(parts))
    function root(i)
        while components[i] != i
            i = components[i]
        end
        return i
    end
    for i in eachindex(parts), j in i+1:length(parts)
        root(i) == root(j) && continue
        connected = false
        for ai in eachindex(parts[i]), bi in eachindex(parts[j])
            a, b = parts[i][ai], parts[i][mod1(ai+1, length(parts[i]))]
            c, d = parts[j][bi], parts[j][mod1(bi+1, length(parts[j]))]
            direction = normalize(b-a)
            if abs(cross2(c-a, direction)) <= tol && abs(cross2(d-a, direction)) <= tol
                lo, hi = minmax(dot(c-a, direction), dot(d-a, direction))
                if min(norm(b-a), hi)-max(0, lo) > tol
                    connected = true
                    break
                end
            end
        end
        if connected
            components[root(j)] = root(i)
        end
    end
    groups = Dict{Int, Vector{Vector{Point2{T}}}}()
    for i in eachindex(parts)
        push!(get!(groups, root(i), Vector{Point2{T}}[]), parts[i])
    end
    return [groups[i] for i in sort!(collect(keys(groups)))]
end

function coincident_faces(mesh, cell, frame, tol)
    out = Int[]
    for face in mesh.faces.cells_to_faces[cell]
        nodes = mesh.faces.faces_to_nodes[face]
        if !isempty(nodes) && all(node -> abs(dot(mesh.node_points[node]-frame.origin, frame.normal)) <= tol, nodes)
            push!(out, face)
        end
    end
    return out
end

function segment_overlap(a, b, c, d, tol)
    L = norm(b-a)
    direction = (b-a)/L
    if norm(cross(c-a, direction)) > tol || norm(cross(d-a, direction)) > tol
        return nothing
    end
    lo, hi = minmax(dot(c-a, direction), dot(d-a, direction))
    lo, hi = max(0, lo), min(L, hi)
    if hi-lo <= tol
        return nothing
    end
    return (a+lo*direction, a+hi*direction)
end

function interface_measure(segments, tol)
    T = typeof(tol)
    # Triangle fans can produce duplicate or partially overlapping segments.
    groups = Vector{Tuple{Point3{T}, Point3{T}}}[]
    for (a, b) in segments
        assigned = false
        for group in groups
            c, d = first(group)
            direction = normalize(d-c)
            if norm(cross(a-c, direction)) <= tol && norm(cross(b-c, direction)) <= tol
                push!(group, (a, b))
                assigned = true
                break
            end
        end
        if !assigned
            push!(groups, [(a, b)])
        end
    end
    length_sum, center = zero(T), zero(Point3{T})
    for group in groups
        origin, endpoint = first(group)
        direction = normalize(endpoint-origin)
        intervals = [minmax(dot(a-origin, direction), dot(b-origin, direction)) for (a, b) in group]
        for (lo, hi) in merge_intervals(intervals, tol)
            L = hi-lo
            length_sum += L
            center += L*(origin+(lo+hi)/2*direction)
        end
    end
    if length_sum > 0
        center /= length_sum
    end
    return length_sum, center
end

"""
    embedded_fracture_mesh(matrix_mesh, fractures; tolerance=nothing, intersections=true, float_type=Jutul.float_type(matrix_mesh))

Clip planar fracture objects against an unmodified 3D mesh. A bounding-box tree
prunes matrix cells and fracture pairs. Warped faces use Jutul's face-center
triangle fans. Fractures on a matrix face share a single fracture cell connected
to both matrix cells. Point contacts have zero transmissibility.

`resolution=(nu, nv)` or `cell_length=h` selects an independent planar fracture
grid. A scalar applies to both axes, a tuple selects each axis, and a vector
supplies one scalar/tuple/`nothing` per input fracture. Overrides replace the
fracture object's `mesh_options`; `nothing` entries retain object settings.
Specify only one option for each fracture. Matrix clipping is the default.
Use `resolution=:matrix` to override an object's independent grid settings.
Grid cells can span multiple matrix cells, using overlap-weighted exchange and
storage. Disconnected active regions are retained as separate fracture cells.
"""
function embedded_fracture_mesh(matrix_mesh::JutulMesh, fractures;
        tolerance = nothing,
        intersections = true,
        float_type = Jutul.float_type(matrix_mesh),
        resolution = nothing,
        cell_length = nothing
    )
    T = float_type
    T <: AbstractFloat || throw(ArgumentError("float_type must be an AbstractFloat type"))
    mesh = UnstructuredMesh(matrix_mesh)
    dim(mesh) == 3 || throw(ArgumentError("EFM requires a three-dimensional matrix mesh"))
    number_of_cells(mesh) > 0 || throw(ArgumentError("EFM requires at least one active matrix cell"))
    fractures = FractureCollection(fractures)
    isempty(fractures) && throw(ArgumentError("at least one fracture is required"))
    index = CellSpatialIndex(mesh, T)
    if isnothing(tolerance)
        root = first(index.nodes)
        scale = max(norm(root.high-root.low), maximum(abs, root.low), maximum(abs, root.high))
        tolerance = max(1e-9, 1e-10*norm(root.high-root.low), 4eps(T)*scale)
    end
    tol = T(tolerance)
    isfinite(tol) && tol > 0 || throw(ArgumentError("tolerance must be positive and finite"))
    frames = fracture_frames(fractures, tol)
    grids = fracture_grids(fractures, frames, tol; resolution = resolution, cell_length = cell_length)
    patches = FracturePatch{T}[]
    face_patches = Dict{Tuple{Int, Int}, Dict{Int, Vector{Tuple{Point3{T}, Point3{T}}}}}()
    coincident = Dict{Tuple{Int, Int}, Vector{Int}}()
    foreach(fractures.groups, fractures.indices) do group, ids
        clip_fracture_group!(patches, face_patches, coincident, mesh, index, group, ids, frames, tol)
    end
    restore_fracture_order!(patches, face_patches)
    patches, face_patches = remesh_fracture_patches(patches, face_patches, grids, frames, mesh, tol)
    restore_fracture_order!(patches, face_patches)
    isempty(patches) && throw(ArgumentError("fractures have no positive-area intersection with the active matrix mesh"))
    pairs = Tuple{Int, Int}[]
    lengths = T[]
    centers = Point3{T}[]
    distances = Tuple{T, T}[]
    types = Symbol[]
    # Shared matrix faces define connections within one fracture. This preserves
    # fault/ACTNUM topology instead of connecting nearby but disconnected cells.
    interfaces = fracture_interfaces(face_patches, tol; mesh = mesh, patches = patches)
    for (pair, segments) in sort!(collect(interfaces); by = first)
        L, center = interface_measure(segments, tol)
        L > tol || continue
        l, r = pair
        frame = frames[patches[l].fracture]
        # TPFA distance normal to the interface, in the fracture plane.
        normal = normalize(cross(last(first(segments))-first(first(segments)), frame.normal))
        dl = abs(dot(patches[l].centroid-center, normal))
        dr = abs(dot(patches[r].centroid-center, normal))
        push!(pairs, pair)
        push!(lengths, L)
        push!(centers, center)
        push!(distances, (max(dl, tol), max(dr, tol)))
        push!(types, :fracture)
    end
    intersection_groups = FractureIntersection{T}[]
    if intersections
        intersection_groups = add_fracture_intersections!(pairs, lengths, centers, distances, types, mesh, patches, frames, interfaces, tol)
    end
    nf, nc = length(pairs), length(patches)
    N = reshape(reinterpret(Int, pairs), 2, nf)
    cc = fracture_point_matrix([p.centroid for p in patches], 3, T)
    fc = fracture_point_matrix(centers, 3, T)
    normals = zeros(T, 3, nf)
    for (i, (l, r)) in enumerate(pairs)
        delta = patches[r].centroid-patches[l].centroid
        if norm(delta) > tol
            normals[:, i] = normalize(delta)
        else
            normals[:, i] = frames[patches[l].fracture].u
        end
    end
    geo = TwoPointFiniteVolumeGeometry(Matrix(N), lengths, [p.area for p in patches], normals, Matrix(cc), Matrix(fc);
        half_face_cells = vec(Matrix(N)), half_face_faces = repeat(collect(1:nf); inner = 2))
    mc, fcells = Int[], Int[]
    mdist, marea = T[], T[]
    matrix_distances = Dict{Tuple{Int, Int}, T}()
    for (fcell, patch) in enumerate(patches)
        frame = frames[patch.fracture]
        for (i, cell) in enumerate(patch.matrix_cells)
            push!(mc, cell)
            push!(fcells, fcell)
            distance = get!(matrix_distances, (patch.fracture, cell)) do
                max(cell_average_distance(mesh, cell, frame.origin, frame.normal), tol)
            end
            push!(mdist, distance)
            push!(marea, patch.matrix_areas[i])
        end
    end
    dist = fracture_point_matrix(distances, 2, T)
    geometry = EmbeddedFractureGeometry(frames, grids, patches, geo, Matrix(dist), types, intersection_groups, mc, fcells, mdist, marea, tol)
    return EmbeddedFractureMesh(matrix_mesh, fractures, geometry, geometry_tags(geo))
end

function fracture_interfaces(face_patches, tol::T; mesh = nothing, patches = nothing) where T
    interfaces = Dict{Tuple{Int, Int}, Vector{Tuple{Point3{T}, Point3{T}}}}()
    for ((fi, face), connected) in pairs(face_patches)
        ids = sort!(collect(keys(connected)))
        for i in 1:length(ids), j in i+1:length(ids)
            l, r = ids[i], ids[j]
            if face < 0 && !isnothing(mesh) &&
                    !connected_matrix_cells(mesh, patches[l].matrix_cells, patches[r].matrix_cells)
                continue
            end
            for (a, b) in connected[l], (c, d) in connected[r]
                segment = segment_overlap(a, b, c, d, tol)
                if !isnothing(segment)
                    push!(get!(interfaces, (l, r), Tuple{Point3{T}, Point3{T}}[]), segment)
                end
            end
        end
    end
    return interfaces
end

function matrix_exchange_area(mesh, cell, frame, area, tol)
    factor = 2
    if !isempty(coincident_faces(mesh, cell, frame, tol))
        factor = 1
    else
        for face in mesh.boundary_faces.cells_to_faces[cell]
            nodes = mesh.boundary_faces.faces_to_nodes[face]
            if !isempty(nodes) && all(n -> abs(dot(mesh.node_points[n]-frame.origin, frame.normal)) <= tol, nodes)
                factor = 1
                break
            end
        end
    end
    return factor*area
end

function connected_matrix_cells(mesh, left, right)
    any(c -> c in right, left) && return true
    for cell in left, face in mesh.faces.cells_to_faces[cell]
        isempty(mesh.faces.faces_to_nodes[face]) && continue
        any(c -> c in right, mesh.faces.neighbors[face]) && return true
    end
    return false
end

function add_fracture_intersections!(pairs, lengths, centers, distances, types, matrix_mesh, patches, frames, interfaces, tol)
    T = typeof(tol)
    index = fracture_patch_spatial_index(patches, frames)
    low, high = index.low, index.high
    # Each group represents one geometric intersection line. Its intervals are
    # split at patch endpoints before eliminating the virtual intersection cell.
    groups = Tuple{Point3{T}, Point3{T}, Vector{Tuple{T, T, Int, Int}}}[]
    fracture_groups = Dict{Tuple{Int, Int}, Int}()
    for (l, pl) in enumerate(patches)
        fl = frames[pl.fracture]
        for r in query_cells(index, low[l], high[l], tol)
            r > l || continue
            pr = patches[r]
            pl.fracture != pr.fracture || continue
            connected_matrix_cells(matrix_mesh, pl.matrix_cells, pr.matrix_cells) || continue
            fr = frames[pr.fracture]
            axis = cross(fl.normal, fr.normal)
            sine = norm(axis)
            if sine < 1e-10
                if abs(dot(fr.origin-fl.origin, fl.normal)) <= tol
                    other = [orient_polygon([project(fl, unproject(fr, p)) for p in part]) for part in pr.parts]
                    overlap = 0
                    for a in pl.parts, b in other
                        p = clip_polygon(a, b, tol)
                        if length(p) >= 3
                            overlap += polygon_area_centroid(p)[1]
                        end
                    end
                    overlap <= tol^2 || throw(ArgumentError("coplanar overlapping fractures $(pl.fracture) and $(pr.fracture) must be merged before preprocessing"))
                end
                continue
            end
            axis /= sine
            tangent = cross(axis, fl.normal)
            origin = fl.origin + dot(fr.origin-fl.origin, fr.normal)/dot(tangent, fr.normal)*tangent
            ol, or = project(fl, origin), project(fr, origin)
            dl, dr = Point2{T}(dot(axis, fl.u), dot(axis, fl.v)), Point2{T}(dot(axis, fr.u), dot(axis, fr.v))
            intervals = overlap_intervals(line_intervals(pl.parts, ol, dl, tol), line_intervals(pr.parts, or, dr, tol), tol)
            isempty(intervals) && continue
            fracture_pair = minmax(pl.fracture, pr.fracture)
            group_id = get(fracture_groups, fracture_pair, 0)
            if group_id == 0
                for (i, (o, d, records)) in enumerate(groups)
                    if norm(cross(d, axis)) <= 1e-10 && norm(cross(origin-o, d)) <= tol
                        group_id = i
                        break
                    end
                end
            end
            if group_id == 0
                push!(groups, (origin, axis, Tuple{T, T, Int, Int}[]))
                group_id = length(groups)
            end
            fracture_groups[fracture_pair] = group_id
            o, d, records = groups[group_id]
            for (a, b) in intervals
                lo, hi = minmax(dot(origin+a*axis-o, d), dot(origin+b*axis-o, d))
                push!(records, (lo, hi, l, r))
            end
        end
    end
    sections = Tuple{Point3{T}, Point3{T}, T, Vector{Int}}[]
    for (origin, axis, records) in groups
        endpoints = sort!(vcat([a for (a, b, l, r) in records], [b for (a, b, l, r) in records]))
        unique_points = T[]
        for p in endpoints
            if isempty(unique_points) || p-last(unique_points) > tol
                push!(unique_points, p)
            end
        end
        for i in 1:length(unique_points)-1
            a, b = unique_points[i], unique_points[i+1]
            b-a > tol || continue
            mid = (a+b)/2
            cells = Int[]
            for (lo, hi, l, r) in records
                if lo-tol <= mid <= hi+tol
                    push!(cells, l, r)
                end
            end
            sort!(unique!(cells))
            length(cells) >= 2 || continue
            push!(sections, (origin+mid*axis, axis, b-a, cells))
        end
    end
    # Replace any ordinary same-fracture connection along an intersection with
    # the corresponding star-delta connection, rather than counting it twice.
    cell_faces = [Int[] for _ in patches]
    for (i, (l, r)) in enumerate(pairs)
        push!(cell_faces[l], i)
        push!(cell_faces[r], i)
    end
    for (center, axis, L, cells) in sections
        a, b = center-L/2*axis, center+L/2*axis
        candidate_faces = unique!(reduce(vcat, cell_faces[cells]))
        for i in candidate_faces
            pair = pairs[i]
            l, r = pair
            if l in cells && r in cells && haskey(interfaces, pair)
                segments = Tuple{Point3{T}, Point3{T}}[]
                for (c, d) in interfaces[pair]
                    overlap = segment_overlap(a, b, c, d, tol)
                    if !isnothing(overlap)
                        push!(segments, overlap)
                    end
                end
                removed, _ = interface_measure(segments, tol)
                lengths[i] = max(lengths[i]-removed, 0)
            end
        end
    end
    keep = findall(>(tol), lengths)
    for data in (pairs, lengths, centers, distances, types)
        kept = data[keep]
        empty!(data)
        append!(data, kept)
    end
    intersections = FractureIntersection{T}[]
    for (center, axis, L, cells) in sections
        ds = T[]
        for cell in cells
            patch = patches[cell]
            frame = frames[patch.fracture]
            origin = project(frame, center)
            direction = Point2{T}(dot(axis, frame.u), dot(axis, frame.v))
            normal = Point2{T}(-direction[2], direction[1])
            avg = polygon_average_distance(patch.parts, origin, normal, tol)
            signed = [dot(p-origin, normal) for part in patch.parts for p in part]
            if minimum(signed) < -tol && maximum(signed) > tol
                avg /= 2 # two sides of an unsplit patch feed the intersection
            end
            push!(ds, max(avg, tol))
        end
        faces = Int[]
        for i in 1:length(cells), j in i+1:length(cells)
            l, r = cells[i], cells[j]
            push!(pairs, (l, r))
            push!(lengths, L)
            push!(centers, center)
            push!(distances, (ds[i], ds[j]))
            push!(types, :intersection)
            push!(faces, length(pairs))
        end
        push!(intersections, FractureIntersection(cells, ds, L, center, faces))
    end
    return intersections
end

function fracture_frames(fractures::FractureCollection, tol::T) where T
    frames = Vector{FractureFrame{T}}(undef, length(fractures))
    foreach(fractures.groups, fractures.indices) do group, ids
        fracture_frames!(frames, group, ids, tol)
    end
    return frames
end

function fracture_frames!(frames, group, ids, tol)
    for (fracture, id) in zip(group, ids)
        frames[id] = fracture_frame(fracture, tol)
    end
    return frames
end

function clip_fracture_group!(patches, face_patches, coincident, mesh, index, group, ids, frames, tol::T) where T
    for (fracture, fi) in zip(group, ids)
        frame = frames[fi]
        vertices = fracture_polygon(fracture)
        low, high = reduce((a, b) -> min.(a, b), vertices), reduce((a, b) -> max.(a, b), vertices)
        for cell in query_cells(index, low, high, tol; origin = frame.origin, normal = frame.normal)
            loops, segments = cell_section(mesh, cell, frame.origin, frame.normal, frame.u, frame.v, tol)
            parts = Vector{Point2{T}}[]
            for loop in loops
                sections = try
                    convex_parts(loop, tol; resolve_winding = true)
                catch err
                    throw(ArgumentError("Could not triangulate fracture $fi in matrix cell $cell: $(sprint(showerror, err)). Section vertices: $loop"))
                end
                for section in sections, boundary in frame.parts
                    p = clip_polygon(section, boundary, tol)
                    if length(p) >= 3 && polygon_area_centroid(p)[1] > tol^2
                        push!(parts, p)
                    end
                end
            end
            isempty(parts) && continue
            cf = coincident_faces(mesh, cell, frame, tol)
            for component in connected_patch_parts(parts, tol)
                area, centroid = patch_geometry(component, frame)
                origin = first(first(component))
                diameter = maximum(norm(p-origin) for part in component for p in part)
                # Area alone does not detect long, numerically zero-width
                # sections along split faces. Apply the geometric length
                # tolerance to their width before creating a storage cell.
                area > tol*diameter || continue
                patch = 0
                for face in cf
                    key = (fi, face)
                    for candidate in get(coincident, key, Int[])
                        p = patches[candidate]
                        if abs(p.area-area) <= tol*sqrt(area) && norm(p.centroid-centroid) <= tol
                            patch = candidate
                            break
                        end
                    end
                    if patch > 0
                        break
                    end
                end
                if patch == 0
                    exchange = matrix_exchange_area(mesh, cell, frame, area, tol)
                    push!(patches, FracturePatch(fi, [cell], component, centroid, area, [exchange], [one(T)]))
                    patch = length(patches)
                    for face in cf
                        push!(get!(coincident, (fi, face), Int[]), patch)
                    end
                else
                    push!(patches[patch].matrix_cells, cell)
                    push!(patches[patch].matrix_areas, matrix_exchange_area(mesh, cell, frame, area, tol))
                    push!(patches[patch].matrix_weights, one(T))
                    fill!(patches[patch].matrix_weights, one(T)/length(patches[patch].matrix_cells))
                end
                for (face, a, b) in clipped_face_segments(frame, segments, tol; parts = component)
                    connected = get!(face_patches, (fi, face), Dict{Int, Vector{Tuple{Point3{T}, Point3{T}}}}())
                    push!(get!(connected, patch, Tuple{Point3{T}, Point3{T}}[]), (a, b))
                end
            end
        end
    end
    return patches
end

function restore_fracture_order!(patches, face_patches)
    issorted(p.fracture for p in patches) && return patches
    order = sortperm(patches; by = p -> p.fracture, alg = Base.Sort.MergeSort)
    renumber = invperm(order)
    patches[:] = patches[order]
    for (key, connected) in pairs(face_patches)
        face_patches[key] = Dict(renumber[cell] => segments for (cell, segments) in pairs(connected))
    end
    return patches
end
