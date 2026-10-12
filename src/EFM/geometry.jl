cross2(a, b) = a[1]*b[2] - a[2]*b[1]

function polygon_area_centroid(p)
    T = eltype(eltype(p))
    origin = first(p)
    area = zero(T)
    centroid = zero(eltype(p))
    for i in 2:length(p)-1
        a = cross2(p[i] - origin, p[i+1] - origin)/2
        area += a
        centroid += a*(origin + p[i] + p[i+1])/3
    end
    if area == 0
        return zero(T), origin
    end
    return abs(area), centroid/area
end

function orient_polygon(p)
    s = sum(cross2(p[i] - p[1], p[i+1] - p[1]) for i in 2:length(p)-1)
    if s < 0
        return reverse(p)
    end
    return p
end

function validate_polygon(p, tol)
    n = length(p)
    for i in 1:n
        j = mod1(i+1, n)
        norm(p[i]-p[j]) > tol || throw(ArgumentError("polygon has repeated adjacent vertices"))
        for k in i+1:n
            l = mod1(k+1, n)
            if j == k || l == i
                continue
            end
            a, b, c, d = p[i], p[j], p[k], p[l]
            s1, s2 = cross2(b-a, c-a), cross2(b-a, d-a)
            s3, s4 = cross2(d-c, a-c), cross2(d-c, b-c)
            overlaps = all(min.(a, b) .<= max.(c, d) .+ tol) && all(min.(c, d) .<= max.(a, b) .+ tol)
            if overlaps && s1*s2 <= 0 && s3*s4 <= 0
                throw(ArgumentError("polygon boundary intersects itself"))
            end
        end
    end
    return nothing
end

function convex_polygon(p, tol)
    n = length(p)
    return all(cross2(p[mod1(i+1, n)]-p[i], p[mod1(i+2, n)]-p[mod1(i+1, n)]) >= -tol^2 for i in 1:n)
end

function clean_polygon(p, tol)
    out = eltype(p)[]
    for pt in p
        if isempty(out) || norm(pt-last(out)) > tol
            push!(out, pt)
        end
    end
    if length(out) > 1 && norm(first(out)-last(out)) <= tol
        pop!(out)
    end
    # Collinear vertices obstruct ear clipping and add no geometry.
    changed = true
    while changed && length(out) > 3
        changed = false
        for i in eachindex(out)
            a = out[i]-out[mod1(i-1, length(out))]
            b = out[mod1(i+1, length(out))]-out[i]
            if abs(cross2(a, b)) <= tol*max(norm(a), norm(b))
                deleteat!(out, i)
                changed = true
                break
            end
        end
    end
    return out
end

function convex_parts(p, tol; resolve_winding = false)
    T = eltype(eltype(p))
    p = orient_polygon(clean_polygon(p, tol))
    if length(p) < 3
        return Vector{Point2{T}}[]
    end
    if resolve_winding
        try
            validate_polygon(p, tol)
        catch err
            err isa ArgumentError || rethrow()
            return polygon_winding_parts(p, tol)
        end
    end
    if convex_polygon(p, tol)
        return [p]
    end
    parts = Vector{Point2{T}}[]
    ix = collect(eachindex(p))
    while length(ix) > 3
        found = false
        for j in eachindex(ix)
            ia, ib, ic = ix[mod1(j-1, length(ix))], ix[j], ix[mod1(j+1, length(ix))]
            a, b, c = p[ia], p[ib], p[ic]
            if cross2(b-a, c-b) <= tol^2
                continue
            end
            contains = false
            for k in ix
                if k == ia || k == ib || k == ic
                    continue
                end
                q = p[k]
                if cross2(b-a, q-a) >= -tol^2 && cross2(c-b, q-b) >= -tol^2 && cross2(a-c, q-c) >= -tol^2
                    contains = true
                    break
                end
            end
            if !contains
                push!(parts, [a, b, c])
                deleteat!(ix, j)
                found = true
                break
            end
        end
        if !found
            if resolve_winding
                return polygon_winding_parts(p, tol)
            end
            error("Could not triangulate planar section; check polygon or mesh geometry")
        end
    end
    push!(parts, p[ix])
    return parts
end

"""
Resolve a section with crossing/retraced face-fan segments into convex regions
with positive winding. Slab endpoints include every segment crossing, so the
trapezoids are exact rather than a convex-hull approximation.
"""
function polygon_winding_parts(p, tol)
    T = eltype(eltype(p))
    events = [q[1] for q in p]
    n = length(p)
    for i in 1:n, j in i+1:n
        a, b = p[i], p[mod1(i+1, n)]
        c, d = p[j], p[mod1(j+1, n)]
        u, v = b-a, d-c
        denominator = cross2(u, v)
        if abs(denominator) <= eps(T)*norm(u)*norm(v)
            continue
        end
        t, s = cross2(c-a, v)/denominator, cross2(c-a, u)/denominator
        if 0 < t < 1 && 0 < s < 1
            push!(events, a[1]+t*u[1])
        end
    end
    sort!(events)
    xs = T[]
    for x in events
        if isempty(xs) || x-last(xs) > tol
            push!(xs, x)
        end
    end
    parts = Vector{Point2{T}}[]
    for k in 1:length(xs)-1
        xl, xr = xs[k], xs[k+1]
        xm = (xl+xr)/2
        edges = Tuple{T, T, T, Int}[]
        for i in 1:n
            a, b = p[i], p[mod1(i+1, n)]
            if min(a[1], b[1]) < xm < max(a[1], b[1])
                slope = (b[2]-a[2])/(b[1]-a[1])
                yl, yr = a[2]+slope*(xl-a[1]), a[2]+slope*(xr-a[1])
                sign = -1
                if b[1] > a[1]
                    sign = 1
                end
                push!(edges, ((yl+yr)/2, yl, yr, sign))
            end
        end
        sort!(edges; by = first)
        winding = 0
        for i in 1:length(edges)-1
            lower, upper = edges[i], edges[i+1]
            winding += lower[4]
            if winding > 0 && upper[1]-lower[1] > tol
                q = clean_polygon([Point2{T}(xl, lower[2]), Point2{T}(xr, lower[3]),
                    Point2{T}(xr, upper[3]), Point2{T}(xl, upper[2])], tol)
                if length(q) >= 3 && polygon_area_centroid(q)[1] > tol^2
                    push!(parts, q)
                end
            end
        end
    end
    return parts
end

function clip_halfplane(p, a, b, tol)
    T = eltype(a)
    out = eltype(p)[]
    isempty(p) && return out
    e = b-a
    l = last(p)
    dl = cross2(e, l-a)
    for r in p
        dr = cross2(e, r-a)
        inside_l, inside_r = dl >= -tol*norm(e), dr >= -tol*norm(e)
        if inside_l != inside_r && abs(dl-dr) > eps(T)
            push!(out, l + clamp(dl/(dl-dr), 0, 1)*(r-l))
        end
        if inside_r
            push!(out, r)
        end
        l, dl = r, dr
    end
    return clean_polygon(out, tol)
end

function clip_polygon(p, boundary, tol)
    T = eltype(eltype(p))
    out = p
    for i in eachindex(boundary)
        out = clip_halfplane(out, boundary[i], boundary[mod1(i+1, length(boundary))], tol)
        if length(out) < 3
            return Point2{T}[]
        end
    end
    return out
end

function point_index!(points, p, tol)
    for (i, q) in enumerate(points)
        if norm(p-q) <= tol
            return i
        end
    end
    push!(points, p)
    return length(points)
end

function triangle_plane_segment(a, b, c, plane::Jutul.PlaneCut, tol)
    T = eltype(plane.point)
    pts = Point3{T}[]
    vertices = (a, b, c)
    distances = map(p -> Jutul.CutCellMeshes.signed_distance(plane, p), vertices)
    for i in 1:3
        j = mod1(i+1, 3)
        p, q = vertices[i], vertices[j]
        dp, dq = distances[i], distances[j]
        if abs(dp) <= tol
            point_index!(pts, p-dp*plane.normal, tol)
        end
        if (dp < -tol && dq > tol) || (dq < -tol && dp > tol)
            point_index!(pts, Jutul.CutCellMeshes.edge_plane_intersection(p, q, plane), tol)
        end
    end
    segment = nothing
    if length(pts) == 2 && norm(pts[1]-pts[2]) > tol
        segment = (pts[1], pts[2])
    end
    return segment
end

"""Intersect a cell boundary using the same face-center triangle fan as Jutul geometry."""
function cell_section(mesh, cell, origin, normal, u, v, tol)
    T = eltype(origin)
    plane = Jutul.PlaneCut(origin, normal)
    points = Point3{T}[]
    edges = Set{Tuple{Int, Int}}()
    face_segments = Tuple{Int, Point3{T}, Point3{T}}[]
    coplanar = Vector{Point2{T}}[]
    for faces in (mesh.faces, mesh.boundary_faces)
        for face in faces.cells_to_faces[cell]
            nodes = faces.faces_to_nodes[face]
            isempty(nodes) && continue # Non-neighboring connection without geometry.
            pts = Point3{T}.(mesh.node_points[nodes])
            center = sum(pts)/length(pts)
            face_id = 0
            if faces === mesh.faces
                face_id = face
            end
            if all(p -> abs(dot(p-origin, normal)) <= tol, pts)
                push!(coplanar, [Point2{T}(dot(p-origin, u), dot(p-origin, v)) for p in pts])
                continue
            end
            for i in eachindex(pts)
                seg = triangle_plane_segment(center, pts[i], pts[mod1(i+1, length(pts))], plane, tol)
                isnothing(seg) && continue
                a, b = seg
                ia, ib = point_index!(points, a, tol), point_index!(points, b, tol)
                push!(edges, minmax(ia, ib))
                if face_id > 0
                    push!(face_segments, (face_id, a, b))
                end
            end
        end
    end
    if !isempty(coplanar)
        return coplanar, face_segments
    end
    loops = Vector{Point2{T}}[]
    while !isempty(edges)
        a, b = pop!(edges)
        ids = [a, b]
        current = b
        while current != a
            next_edge = nothing
            for e in edges
                if current in e
                    next_edge = e
                    break
                end
            end
            if isnothing(next_edge)
                # Tangency at a vertex or an edge has no fracture area.
                empty!(ids)
                break
            end
            delete!(edges, next_edge)
            l, r = next_edge
            if l == current
                current = r
            else
                current = l
            end
            if current != a
                push!(ids, current)
            end
        end
        if length(ids) >= 3
            push!(loops, [Point2{T}(dot(points[i]-origin, u), dot(points[i]-origin, v)) for i in ids])
        end
    end
    return loops, face_segments
end

function merge_intervals(intervals, tol)
    T = typeof(tol)
    isempty(intervals) && return intervals
    sort!(intervals; by = first)
    out = Tuple{T, T}[]
    a, b = first(intervals)
    for (c, d) in intervals[2:end]
        if c <= b + tol
            b = max(b, d)
        else
            push!(out, (a, b))
            a, b = c, d
        end
    end
    push!(out, (a, b))
    return out
end

function line_intervals(parts, origin, direction, tol; limits = (-Inf, Inf))
    T = eltype(origin)
    intervals = Tuple{T, T}[]
    for p in parts
        lo, hi = T.(limits)
        for i in eachindex(p)
            a, b = p[i], p[mod1(i+1, length(p))]
            e = b-a
            offset, slope = cross2(e, origin-a), cross2(e, direction)
            if abs(slope) <= tol*norm(e)
                if offset < -tol*norm(e)
                    hi = -Inf
                    break
                end
            elseif slope > 0
                lo = max(lo, -offset/slope)
            else
                hi = min(hi, -offset/slope)
            end
        end
        if hi-lo > tol
            push!(intervals, (lo, hi))
        end
    end
    return merge_intervals(intervals, tol)
end

function overlap_intervals(a, b, tol)
    T = typeof(tol)
    out = Tuple{T, T}[]
    for (al, ar) in a, (bl, br) in b
        l, r = max(al, bl), min(ar, br)
        if r-l > tol
            push!(out, (l, r))
        end
    end
    return merge_intervals(out, tol)
end

function polygon_average_distance(parts, origin, normal, tol)
    T = eltype(origin)
    integral = zero(T)
    area = zero(T)
    direction = Point2{T}(normal[2], -normal[1])
    for p in parts
        A, c = polygon_area_centroid(p)
        area += A
        total = A*dot(c-origin, normal)
        positive = clip_halfplane(p, origin, origin+direction, tol)
        pos = zero(T)
        if length(positive) >= 3
            Ap, cp = polygon_area_centroid(positive)
            pos = Ap*dot(cp-origin, normal)
        end
        integral += 2pos-total
    end
    return max(integral/area, 0)
end

tetra_volume(a, b, c, d) = dot(b-a, cross(c-a, d-a))/6

"""Exact mean absolute plane distance in a tetrahedron, including a plane through vertices."""
function tetra_average_distance(pts, distances)
    positive = findall(>(0), distances)
    negative = findall(<=(0), distances)
    total = sum(distances)/4
    if isempty(positive) || isempty(negative)
        return abs(total)
    end
    if length(positive) == 3
        return tetra_average_distance(pts, -distances)
    end
    V = abs(tetra_volume(pts...))
    V == 0 && return 0
    crossing(i, j) = pts[i] + distances[i]/(distances[i]-distances[j])*(pts[j]-pts[i])
    if length(positive) == 1
        i = only(positive)
        q = map(j -> crossing(i, j), negative)
        pos = abs(tetra_volume(pts[i], q...))/V*distances[i]/4
    else
        i, j = positive
        k, l = negative
        a, b = pts[i], pts[j]
        ac, ad, bc, bd = crossing(i, k), crossing(i, l), crossing(j, k), crossing(j, l)
        pos = abs(tetra_volume(a, b, ac, ad))*(distances[i]+distances[j])/4
        pos += abs(tetra_volume(b, ac, ad, bd))*distances[j]/4
        pos += abs(tetra_volume(b, ac, bc, bd))*distances[j]/4
        pos /= V
    end
    return max(2pos-total, 0)
end

function cell_average_distance(mesh, cell, origin, normal)
    T = eltype(origin)
    # Signed tetrahedra preserve integrals for non-convex, star-shaped cells.
    center = Point3{T}(Jutul.compute_centroid_and_measure(mesh, Cells(), cell)[1])
    volume = zero(T)
    integral = zero(T)
    for faces in (mesh.faces, mesh.boundary_faces)
        for face in faces.cells_to_faces[cell]
            nodes = faces.faces_to_nodes[face]
            isempty(nodes) && continue
            pts = Point3{T}.(mesh.node_points[nodes])
            fc = sum(pts)/length(pts)
            sign = Jutul.face_orientation_sign(faces, face, cell)
            for i in eachindex(pts)
                tet = (center, fc, pts[i], pts[mod1(i+1, length(pts))])
                V = sign*tetra_volume(tet...)
                distances = SVector{4}(map(p -> dot(p-origin, normal), tet))
                volume += V
                integral += V*tetra_average_distance(tet, distances)
            end
        end
    end
    abs(volume) > 0 || throw(ArgumentError("matrix cell $cell has zero volume"))
    return max(integral/volume, 0)
end
