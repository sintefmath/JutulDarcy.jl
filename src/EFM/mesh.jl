"""
A clipped planar fracture cell with geometric quantities stored in `T`.

Fields:
- `fracture`: One-based input fracture ID, indexing `geometry.frames`.
- `matrix_cells`: Original matrix cell IDs touching the patch. Face-coincident
  patches share one fracture cell between their adjacent matrix cells.
- `parts`: Nonoverlapping convex polygons in the fracture frame's 2D coordinates.
  Disconnected regions are stored as separate patches.
- `centroid`: Area-weighted 3D centroid in metres.
- `area`: Surface area in m² before multiplying by physical aperture.
- `matrix_areas`: Exchange areas (m²) aligned with `matrix_cells`, accounting
  for each cell's overlap and the number of sides feeding the fracture.
- `matrix_weights`: Physical fracture-volume fractions assigned to
  `matrix_cells`. They sum to one; coincident cells share their overlap storage.
"""
struct FracturePatch{T <: AbstractFloat}
    fracture::Int
    matrix_cells::Vector{Int}
    parts::Vector{Vector{Point2{T}}}
    centroid::Point3{T}
    area::T
    matrix_areas::Vector{T}
    matrix_weights::Vector{T}
end

"""
Independent internal grid for one fracture, with coordinates stored in `T`.

Fields:
- `origin`: Lower grid corner in the fracture frame's local 2D coordinates (m).
- `axes`: Orthonormal grid axes in the same local coordinates.
- `spacing`: Grid cell lengths along the two axes (m).
- `resolution`: Cell counts along the axes before clipping to the active mesh.
"""
struct FractureGrid{T <: AbstractFloat}
    origin::Point2{T}
    axes::Tuple{Point2{T}, Point2{T}}
    spacing::Tuple{T, T}
    resolution::Tuple{Int, Int}
end

"""
Alias for [`Jutul.CutCellMeshes.PlanarFrame`](@ref), sharing its fields and
constructors. EFM frame coordinates use metres and `parts` describes the full
input fracture in local 2D coordinates. Existing field access is unchanged.
"""
const FractureFrame = Jutul.CutCellMeshes.PlanarFrame

"""
One positive-length intersection interval eliminated by star–delta coupling.

Fields:
- `cells`: Participating EFM fracture cell IDs, indexing `geometry.patches`.
- `distances`: Effective in-plane distances in metres, aligned with `cells`.
  For an unsplit patch feeding the line from both sides, the average absolute
  distance is halved.
- `length`: Intersection interval length in metres.
- `center`: Interval midpoint in 3D, in metres.
- `faces`: Fracture–fracture face IDs created by eliminating this interval.
"""
struct FractureIntersection{T <: AbstractFloat}
    cells::Vector{Int}
    distances::Vector{T}
    length::T
    center::Point3{T}
    faces::Vector{Int}
end

"""
Clipped surfaces, intersections and connection geometry produced by EFM
preprocessing. All geometric coordinates, distances and measures use `T`.
`finite_volume` contains the lower-dimensional TPFA geometry; physical aperture
is applied by `embedded_fracture_domain`.

Fields:
- `frames`: Local frames in input fracture order.
- `grids`: Effective internal grids in input fracture order; `nothing` entries
  use matrix-clipped patches. Records the actual resolution and spacing.
- `patches`: Clipped cells in fracture-cell order, grouped by input fracture ID.
- `finite_volume`: TPFA connection graph. Its `volumes` are patch areas (m²)
  and `areas` are interface lengths (m), before aperture scaling. Normals are
  centroid-to-centroid directions; explicit coefficients use the distances below.
- `connection_distances`: Two rows of left/right in-plane distances (m), with
  columns in `finite_volume.neighbors` order.
- `connection_types`: Per-face `:fracture` or `:intersection` classification.
- `intersections`: Star–delta records for positive-length crossing intervals.
- `matrix_cells`, `fracture_cells`: Matrix–fracture exchange endpoints, using
  local one-based cell IDs in their respective meshes.
- `matrix_distances`: Volume-averaged absolute matrix-to-plane distances (m)
  in exchange order, bounded below by `tolerance`.
- `matrix_areas`: Exchange areas (m²) in the same order: twice the overlap area
  for an interior overlap, once per adjacent cell for a face-coincident overlap.
- `tolerance`: Absolute geometric tolerance in metres.
"""
struct EmbeddedFractureGeometry{T <: AbstractFloat}
    frames::Vector{FractureFrame{T}}
    grids::Vector{Union{Nothing, FractureGrid{T}}}
    patches::Vector{FracturePatch{T}}
    finite_volume::TwoPointFiniteVolumeGeometry
    connection_distances::Matrix{T}
    connection_types::Vector{Symbol}
    intersections::Vector{FractureIntersection{T}}
    matrix_cells::Vector{Int}
    fracture_cells::Vector{Int}
    matrix_distances::Vector{T}
    matrix_areas::Vector{T}
    tolerance::T
end

"""
The lower-dimensional EFM mesh. Geometry and connection data are stored in
`geometry::EmbeddedFractureGeometry{T}`. The original matrix mesh is retained
by reference and is never modified.

Fields:
- `matrix_mesh`: Original 3D matrix mesh, shared by reference.
- `fractures`: A `FractureCollection` preserving input fracture IDs.
- `geometry`: Clipping, intersection and exchange geometry in type `T`.
- `tags`: Jutul entity tags for the lower-dimensional connection graph.
"""
struct EmbeddedFractureMesh{T <: AbstractFloat, M, F} <: Jutul.FiniteVolumeMesh
    matrix_mesh::M
    fractures::F
    geometry::EmbeddedFractureGeometry{T}
    tags::Jutul.MeshEntityTags
end

"""
A unified connection graph with matrix cells followed by EFM or DFM fracture
cells. All geometric floating-point arrays use `T`.

Fields:
- `matrix_mesh`: Matrix mesh, unchanged for EFM and already cut for DFM.
- `fracture_mesh`: `EmbeddedFractureMesh` for EFM or `Jutul.EmbeddedMesh` for DFM.
- `geometry`: Physical TPFA geometry: cell volumes in m³, face areas in m²,
  coordinates in metres, and matrix boundary geometry. Cells are ordered
  `[matrix; fractures]`; faces are `[matrix; fracture–fracture; matrix–fracture]`.
  Added connections use explicit coefficients rather than geometric TPFA.
- `tags`: Entity tags for the unified graph.

`UnstructuredMesh(mesh)` creates 3D prisms for plotting, preserving cell and
connection indices. The simulation domain retains the physical coefficients.
"""
struct EmbeddedReservoirMesh{T <: AbstractFloat, M, F} <: Jutul.FiniteVolumeMesh
    matrix_mesh::M
    fracture_mesh::F
    geometry::TwoPointFiniteVolumeGeometry
    tags::Jutul.MeshEntityTags
end

function EmbeddedReservoirMesh(matrix, fractures, geometry, tags)
    T = Jutul.float_type(fractures)
    if eltype(geometry.volumes) != T
        geometry = convert_fracture_geometry(T, geometry)
    end
    return EmbeddedReservoirMesh{T, typeof(matrix), typeof(fractures)}(matrix, fractures, geometry, tags)
end

function convert_fracture_geometry(::Type{T}, geometry) where T
    numeric(x::Missing) = missing
    numeric(x::AbstractArray) = T.(x)
    return TwoPointFiniteVolumeGeometry(geometry.neighbors, numeric(geometry.areas),
        numeric(geometry.volumes), numeric(geometry.normals),
        numeric(geometry.cell_centroids), numeric(geometry.face_centroids);
        half_face_cells = geometry.half_face_cells,
        half_face_faces = geometry.half_face_faces,
        boundary_areas = numeric(geometry.boundary_areas),
        boundary_centroids = numeric(geometry.boundary_centroids),
        boundary_normals = numeric(geometry.boundary_normals),
        boundary_neighbors = geometry.boundary_neighbors)
end

function fracture_point_matrix(points, dimension, ::Type{T}) where T
    values = Matrix{T}(undef, dimension, length(points))
    for (i, point) in enumerate(points), j in 1:dimension
        values[j, i] = point[j]
    end
    return values
end

const EFMMesh = Union{EmbeddedFractureMesh, EmbeddedReservoirMesh}
Jutul.dim(::EFMMesh) = 3
Jutul.tpfv_geometry(mesh::EmbeddedFractureMesh) = mesh.geometry.finite_volume
Jutul.tpfv_geometry(mesh::EmbeddedReservoirMesh) = mesh.geometry
Jutul.float_type(::EmbeddedFractureMesh{T}) where T = T
Jutul.float_type(::EmbeddedReservoirMesh{T}) where T = T
Jutul.number_of_cells(mesh::EFMMesh) = length(tpfv_geometry(mesh).volumes)
Jutul.number_of_faces(mesh::EFMMesh) = size(tpfv_geometry(mesh).neighbors, 2)

function Jutul.number_of_boundary_faces(mesh::EFMMesh)
    b = tpfv_geometry(mesh).boundary_neighbors
    if ismissing(b)
        return 0
    end
    return length(b)
end

function Jutul.get_neighborship(mesh::EFMMesh; internal = true)
    if internal
        return tpfv_geometry(mesh).neighbors
    end
    return tpfv_geometry(mesh).boundary_neighbors
end

Jutul.mesh_z_is_depth(mesh::EFMMesh) = Jutul.mesh_z_is_depth(mesh.matrix_mesh)
Jutul.mesh_entity_tags(mesh::EFMMesh) = mesh.tags

function geometry_tags(geo)
    nb = 0
    if !ismissing(geo.boundary_neighbors)
        nb = length(geo.boundary_neighbors)
    end
    return Jutul.MeshEntityTags([
        (entity = Cells(), count = length(geo.volumes)),
        (entity = Faces(), count = size(geo.neighbors, 2)),
        (entity = BoundaryFaces(), count = nb),
        (entity = Jutul.HalfFaces(), count = 2size(geo.neighbors, 2))])
end

function Jutul.cell_dims(mesh::EmbeddedFractureMesh, cell)
    patch = mesh.geometry.patches[cell]
    frame = mesh.geometry.frames[patch.fracture]
    points = [unproject(frame, p) for part in patch.parts for p in part]
    low, high = reduce((a, b) -> min.(a, b), points), reduce((a, b) -> max.(a, b), points)
    # Physical thickness prevents zero dimensions for axis-aligned surfaces.
    a = mesh.fractures[patch.fracture].properties.aperture
    return Tuple(high-low + a*abs.(frame.normal))
end

function Jutul.cell_dims(mesh::EmbeddedReservoirMesh, cell)
    nm = number_of_cells(mesh.matrix_mesh)
    if cell <= nm
        return Jutul.cell_dims(mesh.matrix_mesh, cell)
    end
    return Jutul.cell_dims(mesh.fracture_mesh, cell-nm)
end

function Jutul.cell_ijk(mesh::EmbeddedReservoirMesh, cell::Integer)
    matrix = UnstructuredMesh(mesh.matrix_mesh)
    nm = number_of_cells(mesh.matrix_mesh)
    if matrix.structure isa CartesianIndex
        parent = cell
        if cell > nm
            parent = fracture_parent_cell(mesh.fracture_mesh, mesh.matrix_mesh, cell-nm)
        end
        ijk = Jutul.cell_ijk(matrix, parent)
    else
        ijk = (cell, 1, 1)
    end
    return ijk
end

function Jutul.grid_dims_ijk(mesh::EmbeddedReservoirMesh)
    matrix = UnstructuredMesh(mesh.matrix_mesh)
    if matrix.structure isa CartesianIndex
        dims = Jutul.grid_dims_ijk(matrix)
    else
        dims = (number_of_cells(mesh), 1, 1)
    end
    return dims
end

function Jutul.triangulate_mesh(mesh::EmbeddedFractureMesh{T}; outer = false, flatten = true, kwarg...) where T
    points = Point3{T}[]
    triangles = NTuple{3, Int}[]
    cells = Int[]
    for (cell, patch) in enumerate(mesh.geometry.patches)
        frame = mesh.geometry.frames[patch.fracture]
        for part in patch.parts
            offset = length(points)
            append!(points, [unproject(frame, p) for p in part])
            append!(cells, fill(cell, length(part)))
            for i in 2:length(part)-1
                push!(triangles, (offset+1, offset+i, offset+i+1))
            end
        end
    end
    # Jutul's plotting interface maps each triangle vertex back to its cell.
    if flatten
        expanded = Point3{T}[]
        expanded_cells = Int[]
        for tri in triangles
            append!(expanded, points[collect(tri)])
            append!(expanded_cells, cells[collect(tri)])
        end
        points = expanded
        cells = expanded_cells
        triangles = [(3i-2, 3i-1, 3i) for i in eachindex(triangles)]
        points = Matrix(fracture_point_matrix(points, 3, T)')
        triangles = Matrix(reshape(reinterpret(Int, triangles), 3, :)')
    end
    mapper = (Cells = values -> values[cells], indices = (Cells = cells,))
    return (points = points, triangulation = triangles, mapper = mapper)
end

function Jutul.triangulate_mesh(mesh::EmbeddedReservoirMesh{T}; outer = false, flatten = true, kwarg...) where T
    matrix = Jutul.triangulate_mesh(UnstructuredMesh(mesh.matrix_mesh); outer = outer, flatten = flatten, kwarg...)
    fractures = Jutul.triangulate_mesh(mesh.fracture_mesh; flatten = flatten, kwarg...)
    if flatten
        points = T.(vcat(matrix.points, fractures.points))
        offset = size(matrix.points, 1)
        triangles = vcat(matrix.triangulation, fractures.triangulation .+ offset)
    else
        points = vcat(Point3{T}.(matrix.points), fractures.points)
        offset = length(matrix.points)
        triangles = vcat(Tuple.(matrix.triangulation),
            [Tuple(t .+ offset) for t in fractures.triangulation])
    end
    cells = vcat(matrix.mapper.indices.Cells,
        fractures.mapper.indices.Cells .+ number_of_cells(mesh.matrix_mesh))
    mapper = (Cells = values -> values[cells], indices = (Cells = cells,))
    return (points = points, triangulation = triangles, mapper = mapper)
end

function Jutul.plot_primitives(mesh::EFMMesh, plot_type; kwarg...)
    if plot_type == :mesh
        return Jutul.triangulate_mesh(mesh; kwarg...)
    end
    return nothing
end
