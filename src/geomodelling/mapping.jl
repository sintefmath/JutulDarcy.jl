"""
    map_property_to_cells(mesh, source::CartesianMesh, values; mapping = :arithmetic, outside = :nearest, weights = missing, info_level = 0, name = missing, extra_out = false)
    map_property_to_cells(domain::DataDomain, source::CartesianMesh, values; kwarg...)

Map a cell-wise property defined on a Cartesian `source` mesh (for example a
geostatistical realization generated on a box grid) onto the cells of `mesh`.
The target mesh can have arbitrary cell shapes and sizes (corner-point,
tensor-product or unstructured meshes), and can be both finer and coarser than
the source mesh. `values` is either a vector with one value per source cell or
a matrix with one column per source cell (e.g. a 3 by n permeability tensor),
in which case each row is mapped separately.

Each target cell is split into tetrahedra spanned by the cell centroid, a face
centroid and two consecutive face nodes (the same decomposition used to compute
cell volumes). The tetrahedra are then cut along the planes of the source mesh
until each piece lies inside a single source cell. This gives the exact overlap
volume between each target cell and the source cells, and the result is a
volume-weighted average of the source values over each target cell. The
computation is multithreaded over the cells of `mesh`.

# Keyword arguments
- `mapping = :arithmetic`: Volume-weighted average used inside each target
  cell. One of `:arithmetic`, `:geometric` or `:harmonic`. Porosity should use
  `:arithmetic`, while `:geometric` is a reasonable default for permeability.
- `outside = :nearest`: How parts of target cells outside the source mesh are
  treated. `:nearest` uses the nearest source cell, `:ignore` disregards those
  parts (cells entirely outside the source get `NaN`) and `:error` throws.
- `weights = missing`: Overlap volumes from a previous call with
  `extra_out = true`. Pass these to map several properties without
  recomputing the geometric overlap.
- `info_level = 0`: Print a mapping summary if positive.
- `name = missing`: Property name used in the summary.
- `extra_out = false`: Also return the sparse matrix of overlap volumes (target
  cells by source cells), i.e. `(mapped_values, weights)`.
"""
function map_property_to_cells(mesh::JutulMesh, source::CartesianMesh, values;
        mapping::Symbol = :arithmetic,
        outside::Symbol = :nearest,
        weights = missing,
        info_level = 0,
        name = missing,
        extra_out = false
    )
    nsource = number_of_cells(source)
    nvalues = size(values, ndims(values))
    nvalues == nsource || throw(ArgumentError("Number of values $nvalues does not match the number of source cells $nsource."))
    if ismissing(weights)
        weights, outside_volume = cartesian_overlap_weights(mesh, source, outside = outside)
    else
        outside_volume = missing
    end
    mapped_values = apply_cell_weights(weights, values, mapping)

    if info_level > 0
        nc = size(weights, 1)
        str = ismissing(name) ? "" : "$name: "
        msg = "$(str)Mapped $nsource source cells to $nc cells \
            (approx. $(round(nsource/max(nc, 1), sigdigits = 3)) source cells per cell) \
            using $mapping averaging."
        if !ismissing(outside_volume)
            total_volume = sum(weights) + (outside === :ignore ? outside_volume : 0.0)
            frac = outside_volume/max(total_volume, eps())
            msg *= " $(round(100*frac, sigdigits = 3))% of the mesh volume is outside the source mesh (outside = :$outside)."
        end
        jutul_message("Mapping", msg)
    end
    if extra_out
        return (mapped_values, weights)
    else
        return mapped_values
    end
end

"""
    map_property_to_cells(mesh, points, values; volumes = nothing, mapping = :arithmetic, cell_lookup = missing, info_level = 0, name = missing, extra_out = false)

Map a property defined on scattered `points` to the cells of `mesh`. `points`
is either a matrix with three rows and one column per point, or a vector of
points. Each point is assigned to the cell that contains it and the values in
each cell are averaged (weighted by `volumes` if given). Cells that do not
contain any point get the value of the point nearest to the cell centroid, so
coarse point sets can also be mapped to finer meshes. Points outside the mesh
are ignored.

If the property is defined on a Cartesian grid, the method taking a
`CartesianMesh` source is more accurate since it accounts for the overlap
between source and target cells.

With `extra_out = true`, `(mapped_values, cell_lookup)` is returned, where
`cell_lookup` gives the enclosing cell of each point (zero if outside). The
lookup can be passed back through `cell_lookup` to map more properties defined
on the same points.
"""
function map_property_to_cells(mesh::JutulMesh, points::Union{AbstractMatrix, AbstractVector}, values;
        volumes = nothing,
        mapping::Symbol = :arithmetic,
        cell_lookup = missing,
        info_level = 0,
        name = missing,
        extra_out = false
    )
    pts = points_to_svectors(points)
    npts = length(pts)
    nvalues = size(values, ndims(values))
    nvalues == npts || throw(ArgumentError("Number of values $nvalues does not match the number of points $npts."))
    if !isnothing(volumes)
        length(volumes) == npts || throw(ArgumentError("Number of volumes $(length(volumes)) does not match the number of points $npts."))
        all(v -> v > 0, volumes) || throw(ArgumentError("Volumes must be strictly positive."))
    end
    mesh_u = UnstructuredMesh(mesh)
    dim(mesh_u) == 3 || throw(ArgumentError("Expected a 3D mesh."))
    geometry = tpfv_geometry(mesh_u)
    if ismissing(cell_lookup)
        cell_lookup = create_cell_lookup(mesh_u, geometry, pts)
    end
    length(cell_lookup) == npts || throw(ArgumentError("cell_lookup must have one entry per point."))

    nc = number_of_cells(mesh_u)
    I = Int[]
    J = Int[]
    V = Float64[]
    covered = falses(nc)
    for (p, cell) in enumerate(cell_lookup)
        cell == 0 && continue
        covered[cell] = true
        push!(I, cell)
        push!(J, p)
        push!(V, isnothing(volumes) ? 1.0 : Float64(volumes[p]))
    end
    nuncovered = nc - count(covered)
    if nuncovered > 0
        uncovered = findall(!, covered)
        centroids = [SVector{3, Float64}(geometry.cell_centroids[:, c]) for c in uncovered]
        nearest = nearest_point_indices(pts, centroids)
        append!(I, uncovered)
        append!(J, nearest)
        append!(V, ones(nuncovered))
    end
    weights = sparse(I, J, V, nc, npts)
    mapped_values = apply_cell_weights(weights, values, mapping)

    if info_level > 0
        ignored = count(isequal(0), cell_lookup)
        str = ismissing(name) ? "" : "$name: "
        jutul_message("Mapping", "$(str)Mapped $npts points to $nc cells \
            using $mapping averaging. $ignored points were outside the mesh and \
            $nuncovered cells contained no points and were assigned the nearest point value.")
    end
    if extra_out
        return (mapped_values, cell_lookup)
    else
        return mapped_values
    end
end

function map_property_to_cells(domain::DataDomain, source, values; kwarg...)
    return map_property_to_cells(physical_representation(domain), source, values; kwarg...)
end

"""
    map_property_to_cells!(domain, source, values, name; kwarg...)
    map_property_to_cells!(domain, source, properties::AbstractDict; mapping = default_property_mapping, kwarg...)

Mutating variant of [`map_property_to_cells`](@ref) that stores the mapped
values in `domain` under `name`. When a dictionary of properties is given, all
entries are mapped using the same geometric overlap, which is only computed
once. `mapping` can then be a `Symbol` used for all properties, or a function
or dictionary that gives the mapping from the property name. By default,
`:permeability` uses `:geometric` and all other properties `:arithmetic`
averaging. Returns the domain.
"""
function map_property_to_cells!(domain::DataDomain, source, values, name::Symbol; kwarg...)
    domain[name] = map_property_to_cells(domain, source, values; name = name, kwarg...)
    return domain
end

function map_property_to_cells!(domain::DataDomain, source, properties::AbstractDict;
        mapping = default_property_mapping,
        kwarg...
    )
    reuse = missing
    for (name, values) in pairs(properties)
        name = Symbol(name)
        m = mapping isa Symbol ? mapping : get_property_mapping(mapping, name)
        if source isa CartesianMesh
            values, reuse = map_property_to_cells(domain, source, values;
                mapping = m, weights = reuse, name = name, extra_out = true, kwarg...)
        else
            values, reuse = map_property_to_cells(domain, source, values;
                mapping = m, cell_lookup = reuse, name = name, extra_out = true, kwarg...)
        end
        domain[name] = values
    end
    return domain
end

function default_property_mapping(name::Symbol)
    if name == :permeability
        return :geometric
    else
        return :arithmetic
    end
end

get_property_mapping(mapping::AbstractDict, name) = get(mapping, name, default_property_mapping(name))
get_property_mapping(mapping, name) = mapping(name)

function apply_cell_weights(weights::AbstractSparseMatrix, values::AbstractVector, mapping::Symbol)
    wsum = vec(sum(weights, dims = 2))
    v = Float64.(values)
    if mapping === :arithmetic
        mapped = (weights*v)./wsum
    elseif mapping === :geometric
        all(x -> x > 0, v) || throw(ArgumentError("Values must be strictly positive for geometric averaging."))
        mapped = exp.((weights*log.(v))./wsum)
    elseif mapping === :harmonic
        all(x -> x > 0, v) || throw(ArgumentError("Values must be strictly positive for harmonic averaging."))
        mapped = wsum./(weights*inv.(v))
    else
        throw(ArgumentError("mapping must be :arithmetic, :geometric or :harmonic, was :$mapping."))
    end
    # Cells without any overlap with the source are undefined
    for (i, w) in enumerate(wsum)
        if w == 0
            mapped[i] = NaN
        end
    end
    return mapped
end

function apply_cell_weights(weights::AbstractSparseMatrix, values::AbstractMatrix, mapping::Symbol)
    out = zeros(size(values, 1), size(weights, 1))
    for i in axes(values, 1)
        out[i, :] = apply_cell_weights(weights, values[i, :], mapping)
    end
    return out
end

"""
    cartesian_overlap_weights(mesh, source::CartesianMesh; outside = :nearest)

Compute a sparse matrix `W` where `W[c, s]` is the volume of the
overlap between cell `c` of `mesh` and cell `s` of `source`. Also returns the
volume of the mesh that lies outside the source mesh.
"""
function cartesian_overlap_weights(mesh::JutulMesh, source::CartesianMesh; outside::Symbol = :nearest)
    outside in (:nearest, :ignore, :error) || throw(ArgumentError("outside must be :nearest, :ignore or :error, was :$outside."))
    mesh_u = UnstructuredMesh(mesh)
    dim(mesh_u) == 3 || throw(ArgumentError("Expected a 3D mesh."))
    length(source.dims) == 3 || throw(ArgumentError("Expected a 3D source mesh."))
    lookup = CartesianLookup(cartesian_node_coordinates(source), source.dims, outside)

    geometry = tpfv_geometry(mesh_u)
    T = SVector{3, Float64}
    cell_centroids = vec(reinterpret(T, geometry.cell_centroids))
    face_centroids = vec(reinterpret(T, geometry.face_centroids))
    boundary_centroids = vec(reinterpret(T, geometry.boundary_centroids))
    nc = number_of_cells(mesh_u)
    # Cells are independent, so chunks of cells are processed in parallel
    nchunks = min(nc, 4*Threads.nthreads())
    chunks = Iterators.partition(1:nc, cld(nc, max(nchunks, 1)))
    tasks = map(chunks) do cells
        Threads.@spawn cartesian_overlap_weights_for_cells(mesh_u, lookup, cells,
            cell_centroids, face_centroids, boundary_centroids, geometry.volumes)
    end
    results = try
        fetch.(tasks)
    catch e
        if e isa TaskFailedException
            throw(e.task.exception)
        else
            rethrow()
        end
    end
    I = reduce(vcat, map(r -> r[1], results))
    J = reduce(vcat, map(r -> r[2], results))
    V = reduce(vcat, map(r -> r[3], results))
    outside_volume = sum(r -> r[4], results)
    weights = sparse(I, J, V, nc, number_of_cells(source))
    return (weights, outside_volume)
end

function cartesian_overlap_weights_for_cells(mesh_u, lookup, cells, cell_centroids, face_centroids, boundary_centroids, volumes)
    pts = mesh_u.node_points
    I = Int[]
    J = Int[]
    V = Float64[]
    cell_weights = Dict{Int, Float64}()
    outside_volume = Ref(0.0)
    for c in cells
        empty!(cell_weights)
        xc = cell_centroids[c]
        lo, hi = cell_bounding_box(mesh_u, c)
        if all(axis -> isnothing(splitting_plane(lookup, lo[axis], hi[axis], axis)), 1:3)
            # The whole cell is inside a single source cell
            add_volume_at_point!(cell_weights, outside_volume, lookup, xc, volumes[c])
        else
            for f in mesh_u.faces.cells_to_faces[c]
                face_nodes = mesh_u.faces.faces_to_nodes[f]
                add_face_tetrahedra!(cell_weights, outside_volume, lookup, xc, face_centroids[f], face_nodes, pts)
            end
            for bf in mesh_u.boundary_faces.cells_to_faces[c]
                face_nodes = mesh_u.boundary_faces.faces_to_nodes[bf]
                add_face_tetrahedra!(cell_weights, outside_volume, lookup, xc, boundary_centroids[bf], face_nodes, pts)
            end
        end
        for (s, w) in cell_weights
            push!(I, c)
            push!(J, s)
            push!(V, w)
        end
    end
    return (I, J, V, outside_volume[])
end

struct CartesianLookup
    nodes::NTuple{3, Vector{Float64}}
    dims::NTuple{3, Int}
    tol::NTuple{3, Float64}
    outside::Symbol
end

function CartesianLookup(nodes, dims, outside)
    tol = map(x -> 1e-10*(x[end] - x[1]), nodes)
    return CartesianLookup(nodes, Tuple(dims), tol, outside)
end

function cartesian_node_coordinates(g::CartesianMesh)
    function axis_nodes(i)
        Δ = g.deltas[i]
        x0 = Float64(g.origin[i])
        n = g.dims[i]
        if Δ isa Real
            return collect(x0 .+ Δ .* (0:n))
        else
            return x0 .+ [0.0; cumsum(Float64.(Δ))]
        end
    end
    return (axis_nodes(1), axis_nodes(2), axis_nodes(3))
end

function add_face_tetrahedra!(cell_weights, outside_volume, lookup, xc, xf, face_nodes, pts)
    n = length(face_nodes)
    for i in 1:n
        a = pts[face_nodes[i]]
        b = pts[face_nodes[mod1(i + 1, n)]]
        add_tetrahedron!(cell_weights, outside_volume, lookup, xc, xf, a, b)
    end
end

"""
Find a source node plane along `axis` strictly between `lo` and `hi`. Returns
the plane coordinate, or `nothing` if the interval is inside a single slab.
"""
function splitting_plane(lookup::CartesianLookup, lo, hi, axis)
    nd = lookup.nodes[axis]
    tol = lookup.tol[axis]
    first_plane = searchsortedlast(nd, lo + tol) + 1
    last_plane = searchsortedfirst(nd, hi - tol) - 1
    if first_plane > last_plane
        return nothing
    else
        return nd[(first_plane + last_plane) ÷ 2]
    end
end

"""
Add the overlap volumes between tetrahedron `(a, b, c, d)` and the source
cells to `cell_weights`. The tetrahedron is split recursively along the source
node planes it crosses until each piece lies inside a single source cell, so
the overlap volumes are exact.
"""
function add_tetrahedron!(cell_weights, outside_volume, lookup, a, b, c, d)
    vol = abs(dot(b - a, cross(c - a, d - a)))/6
    vol > 0 || return
    lo = min.(a, b, c, d)
    hi = max.(a, b, c, d)
    for axis in 1:3
        p = splitting_plane(lookup, lo[axis], hi[axis], axis)
        if !isnothing(p)
            split_tetrahedron!(cell_weights, outside_volume, lookup, a, b, c, d, axis, p)
            return
        end
    end
    # Inside a single source cell, or outside the source mesh
    add_volume_at_point!(cell_weights, outside_volume, lookup, (a + b + c + d)/4, vol)
    return
end

"""
Add volume `vol` to the source cell containing `x`, which is assumed to be
representative for the whole volume.
"""
function add_volume_at_point!(cell_weights, outside_volume, lookup, x, vol)
    nx, ny, nz = lookup.dims
    xn, yn, zn = lookup.nodes
    i = clamp(searchsortedlast(xn, x[1]), 1, nx)
    j = clamp(searchsortedlast(yn, x[2]), 1, ny)
    k = clamp(searchsortedlast(zn, x[3]), 1, nz)
    inside = xn[1] <= x[1] <= xn[end] && yn[1] <= x[2] <= yn[end] && zn[1] <= x[3] <= zn[end]
    if !inside
        if lookup.outside === :error
            throw(ArgumentError("Mesh extends outside the source mesh. Use outside = :nearest or :ignore to allow this."))
        end
        outside_volume[] += vol
        lookup.outside === :ignore && return
    end
    idx = (k - 1)*nx*ny + (j - 1)*nx + i
    cell_weights[idx] = get(cell_weights, idx, 0.0) + vol
    return
end

function split_tetrahedron!(cell_weights, outside_volume, lookup, a, b, c, d, axis, p)
    tol = lookup.tol[axis]
    v = SVector(a, b, c, d)
    side = map(x -> x[axis] > p + tol ? 1 : (x[axis] < p - tol ? -1 : 0), v)
    # Order vertices as positive, on the plane, negative
    perm = sortperm(side, rev = true)
    v = v[perm]
    np = count(isequal(1), side)
    nn = count(isequal(-1), side)
    function cut(u, w)
        x = u + ((p - u[axis])/(w[axis] - u[axis]))*(w - u)
        return Base.setindex(x, p, axis)
    end
    add!(t1, t2, t3, t4) = add_tetrahedron!(cell_weights, outside_volume, lookup, t1, t2, t3, t4)
    function add_prism!(x1, x2, x3, y1, y2, y3)
        add!(x1, x2, x3, y3)
        add!(x1, x2, y3, y2)
        add!(x1, y2, y3, y1)
    end
    if np == 1 && nn == 1
        pos, z1, z2, neg = v
        cm = cut(pos, neg)
        add!(pos, z1, z2, cm)
        add!(neg, z1, z2, cm)
    elseif np == 1 && nn == 2
        pos, z, n1, n2 = v
        c1, c2 = cut(pos, n1), cut(pos, n2)
        add!(pos, c1, c2, z)
        add!(z, n1, n2, c2)
        add!(z, n1, c2, c1)
    elseif np == 2 && nn == 1
        p1, p2, z, neg = v
        c1, c2 = cut(p1, neg), cut(p2, neg)
        add!(neg, c1, c2, z)
        add!(z, p1, p2, c2)
        add!(z, p1, c2, c1)
    elseif np == 1 && nn == 3
        pos, n1, n2, n3 = v
        c1, c2, c3 = cut(pos, n1), cut(pos, n2), cut(pos, n3)
        add!(pos, c1, c2, c3)
        add_prism!(n1, n2, n3, c1, c2, c3)
    elseif np == 3 && nn == 1
        p1, p2, p3, neg = v
        c1, c2, c3 = cut(p1, neg), cut(p2, neg), cut(p3, neg)
        add!(neg, c1, c2, c3)
        add_prism!(p1, p2, p3, c1, c2, c3)
    elseif np == 2 && nn == 2
        p1, p2, n1, n2 = v
        c11, c12 = cut(p1, n1), cut(p1, n2)
        c21, c22 = cut(p2, n1), cut(p2, n2)
        add_prism!(p1, c11, c12, p2, c21, c22)
        add_prism!(n1, c11, c21, n2, c12, c22)
    else
        error("Unexpected split of tetrahedron with $np vertices above and $nn below plane.")
    end
    return
end

function points_to_svectors(points::AbstractMatrix)
    size(points, 1) == 3 || throw(ArgumentError("Point matrix must have 3 rows (x, y, z), got $(size(points, 1))."))
    return [SVector{3, Float64}(points[1, j], points[2, j], points[3, j]) for j in axes(points, 2)]
end

function points_to_svectors(points::AbstractVector)
    out = Vector{SVector{3, Float64}}(undef, length(points))
    for (i, pt) in enumerate(points)
        length(pt) == 3 || throw(ArgumentError("Points must have three coordinates (x, y, z), got $(length(pt))."))
        out[i] = SVector{3, Float64}(pt[1], pt[2], pt[3])
    end
    return out
end

"""
Uniform grid of buckets over a bounding box, used to find candidate cells or
points near a given point.
"""
struct BucketGrid
    low::SVector{3, Float64}
    h::SVector{3, Float64}
    dims::NTuple{3, Int}
    buckets::Vector{Vector{Int}}
end

function BucketGrid(low, high, n_items)
    nb = max(1, ceil(Int, cbrt(n_items)))
    h = max.(high - low, 1e-12)./nb
    dims = (nb, nb, nb)
    buckets = [Int[] for _ in 1:prod(dims)]
    return BucketGrid(low, h, dims, buckets)
end

function bucket_ijk(bg::BucketGrid, x)
    return ntuple(d -> clamp(floor(Int, (x[d] - bg.low[d])/bg.h[d]) + 1, 1, bg.dims[d]), 3)
end

function bucket_index(bg::BucketGrid, ijk)
    nx, ny, _ = bg.dims
    return (ijk[3] - 1)*nx*ny + (ijk[2] - 1)*nx + ijk[1]
end

function cell_bounding_box(mesh::UnstructuredMesh, c)
    T = SVector{3, Float64}
    lo = T(Inf, Inf, Inf)
    hi = T(-Inf, -Inf, -Inf)
    for faces in (mesh.faces, mesh.boundary_faces)
        for f in faces.cells_to_faces[c]
            for node in faces.faces_to_nodes[f]
                x = mesh.node_points[node]
                lo = min.(lo, x)
                hi = max.(hi, x)
            end
        end
    end
    return (lo, hi)
end

function create_cell_lookup(mesh::UnstructuredMesh, geometry, pts)
    T = SVector{3, Float64}
    normals = vec(reinterpret(T, geometry.normals))
    face_centroids = vec(reinterpret(T, geometry.face_centroids))
    boundary_centroids = vec(reinterpret(T, geometry.boundary_centroids))
    boundary_normals = vec(reinterpret(T, geometry.boundary_normals))

    # Bin the cells into buckets using their bounding boxes so that only a few
    # candidate cells must be checked for each point.
    nc = number_of_cells(mesh)
    boxes = map(c -> cell_bounding_box(mesh, c), 1:nc)
    low = reduce((a, b) -> min.(a, b), first.(boxes))
    high = reduce((a, b) -> max.(a, b), last.(boxes))
    bg = BucketGrid(low, high, nc)
    for (c, (lo, hi)) in enumerate(boxes)
        blo = bucket_ijk(bg, lo)
        bhi = bucket_ijk(bg, hi)
        for k in blo[3]:bhi[3], j in blo[2]:bhi[2], i in blo[1]:bhi[1]
            push!(bg.buckets[bucket_index(bg, (i, j, k))], c)
        end
    end
    cell_lookup = zeros(Int, length(pts))
    for (idx, pt) in enumerate(pts)
        if any(pt .< low) || any(pt .> high)
            continue
        end
        candidates = bg.buckets[bucket_index(bg, bucket_ijk(bg, pt))]
        cell = Jutul.find_enclosing_cell(mesh, pt, normals, face_centroids, boundary_normals, boundary_centroids, candidates)
        if !isnothing(cell)
            cell_lookup[idx] = cell
        end
    end
    return cell_lookup
end

"""
Index of the nearest point in `pts` for each point in `queries`.
"""
function nearest_point_indices(pts::Vector{SVector{3, Float64}}, queries::Vector{SVector{3, Float64}})
    low = reduce((a, b) -> min.(a, b), pts)
    high = reduce((a, b) -> max.(a, b), pts)
    bg = BucketGrid(low, high, length(pts))
    for (i, pt) in enumerate(pts)
        push!(bg.buckets[bucket_index(bg, bucket_ijk(bg, pt))], i)
    end
    hmin = minimum(bg.h)
    rmax = maximum(bg.dims)
    nearest = zeros(Int, length(queries))
    for (q, x) in enumerate(queries)
        center = bucket_ijk(bg, x)
        best = 0
        best_dist = Inf
        r = 0
        while true
            # Search the shell of buckets at Chebyshev distance r
            for k in (center[3] - r):(center[3] + r), j in (center[2] - r):(center[2] + r), i in (center[1] - r):(center[1] + r)
                if max(abs(i - center[1]), abs(j - center[2]), abs(k - center[3])) != r
                    continue
                end
                if !(1 <= i <= bg.dims[1] && 1 <= j <= bg.dims[2] && 1 <= k <= bg.dims[3])
                    continue
                end
                for p in bg.buckets[bucket_index(bg, (i, j, k))]
                    d = norm(pts[p] - x)
                    if d < best_dist
                        best_dist = d
                        best = p
                    end
                end
            end
            # Points in the remaining shells are at least r*hmin away
            if (best > 0 && best_dist <= r*hmin) || r >= rmax
                break
            end
            r += 1
        end
        nearest[q] = best
    end
    return nearest
end

"""
    mesh_has_ijk(mesh)

Check if the cells of a 3D mesh have a logical IJK index.
"""
mesh_has_ijk(mesh) = false
mesh_has_ijk(mesh::CartesianMesh) = length(mesh.dims) == 3
mesh_has_ijk(mesh::UnstructuredMesh{3, CartesianIndex{3}}) = true

"""
    cell_ijk_linear_indices(mesh)

Linear index `(k-1)*nx*ny + (j-1)*nx + i` of each cell in the logical IJK grid
of the mesh.
"""
function cell_ijk_linear_indices(mesh)
    nx, ny, nz = grid_dims_ijk(mesh)
    return map(1:number_of_cells(mesh)) do c
        i, j, k = cell_ijk(mesh, c)
        (k - 1)*nx*ny + (j - 1)*nx + i
    end
end

"""
    average_cell_spacing(mesh)

Approximate the average cell size in each direction. The horizontal sizes are
the mean extents of the cell bounding boxes, while the vertical size is the
mean of the cell volume divided by the horizontal extents, which is less
sensitive to dipping layers than the vertical extent of the bounding box.
"""
function average_cell_spacing(mesh)
    mesh_u = UnstructuredMesh(mesh)
    geometry = tpfv_geometry(mesh_u)
    nc = number_of_cells(mesh_u)
    dx = dy = dz = 0.0
    for c in 1:nc
        lo, hi = cell_bounding_box(mesh_u, c)
        ex, ey = hi[1] - lo[1], hi[2] - lo[2]
        dx += ex
        dy += ey
        dz += geometry.volumes[c]/max(ex*ey, eps())
    end
    return (dx/nc, dy/nc, dz/nc)
end

"""
    mesh_bounding_box(mesh)

Lower and upper corner of the bounding box of the mesh nodes.
"""
function mesh_bounding_box(mesh)
    pts = UnstructuredMesh(mesh).node_points
    low = reduce((a, b) -> min.(a, b), pts)
    high = reduce((a, b) -> max.(a, b), pts)
    return (low, high)
end
