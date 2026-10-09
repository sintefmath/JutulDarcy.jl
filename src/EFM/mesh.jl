"""A clipped planar fracture cell; `parts` are nonoverlapping convex polygons."""
struct FracturePatch
    fracture::Int
    matrix_cells::Vector{Int}
    parts::Vector{Vector{Point2}}
    centroid::Point3
    area::Float64
end

struct FractureFrame
    origin::Point3
    normal::Point3
    u::Point3
    v::Point3
    parts::Vector{Vector{Point2}}
end

struct FractureIntersection
    cells::Vector{Int}
    distances::Vector{Float64}
    length::Float64
    center::Point3
    faces::Vector{Int}
end

project(frame, p) = Point2(dot(p-frame.origin, frame.u), dot(p-frame.origin, frame.v))
unproject(frame, p) = frame.origin + p[1]*frame.u + p[2]*frame.v

"""
The lower-dimensional EFM mesh. Its cells are fracture polygons clipped to
matrix cells. Faces represent fracture connections, including intersections.
The matrix mesh is retained by reference and is never modified.
"""
struct EmbeddedFractureMesh{M} <: Jutul.FiniteVolumeMesh
    matrix_mesh::M
    fractures::Vector{AbstractEmbeddedFracture}
    frames::Vector{FractureFrame}
    patches::Vector{FracturePatch}
    geometry::TwoPointFiniteVolumeGeometry
    connection_distances::Matrix{Float64}
    connection_types::Vector{Symbol}
    intersections::Vector{FractureIntersection}
    matrix_cells::Vector{Int}
    fracture_cells::Vector{Int}
    matrix_distances::Vector{Float64}
    matrix_areas::Vector{Float64}
    tolerance::Float64
    tags::Jutul.MeshEntityTags
end

"""A connection graph with the original matrix cells followed by EFM cells."""
struct EmbeddedReservoirMesh{M, F} <: Jutul.FiniteVolumeMesh
    matrix_mesh::M
    fracture_mesh::F
    geometry::TwoPointFiniteVolumeGeometry
    tags::Jutul.MeshEntityTags
end

const EFMMesh = Union{EmbeddedFractureMesh, EmbeddedReservoirMesh}
Jutul.dim(::EFMMesh) = 3
Jutul.tpfv_geometry(mesh::EFMMesh) = mesh.geometry
Jutul.number_of_cells(mesh::EFMMesh) = length(mesh.geometry.volumes)
Jutul.number_of_faces(mesh::EFMMesh) = size(mesh.geometry.neighbors, 2)

function Jutul.number_of_boundary_faces(mesh::EFMMesh)
    b = mesh.geometry.boundary_neighbors
    if ismissing(b)
        return 0
    end
    return length(b)
end

function Jutul.get_neighborship(mesh::EFMMesh; internal = true)
    if internal
        return mesh.geometry.neighbors
    end
    return mesh.geometry.boundary_neighbors
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
    patch = mesh.patches[cell]
    frame = mesh.frames[patch.fracture]
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
    nm = number_of_cells(mesh.matrix_mesh)
    if cell <= nm
        return Jutul.cell_ijk(mesh.matrix_mesh, cell)
    end
    return Jutul.cell_ijk(mesh.matrix_mesh, first(mesh.fracture_mesh.patches[cell-nm].matrix_cells))
end

Jutul.grid_dims_ijk(mesh::EmbeddedReservoirMesh) = Jutul.grid_dims_ijk(mesh.matrix_mesh)

function Jutul.triangulate_mesh(mesh::EmbeddedFractureMesh; outer = false, flatten = true, kwarg...)
    points = Point3[]
    triangles = NTuple{3, Int}[]
    cells = Int[]
    for (cell, patch) in enumerate(mesh.patches)
        frame = mesh.frames[patch.fracture]
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
        expanded = Point3[]
        expanded_cells = Int[]
        for tri in triangles
            append!(expanded, points[collect(tri)])
            append!(expanded_cells, cells[collect(tri)])
        end
        points = expanded
        cells = expanded_cells
        triangles = [(3i-2, 3i-1, 3i) for i in eachindex(triangles)]
        points = Matrix(reshape(reinterpret(Float64, points), 3, :)')
        triangles = Matrix(reshape(reinterpret(Int, triangles), 3, :)')
    end
    mapper = (Cells = values -> values[cells], indices = (Cells = cells,))
    return (points = points, triangulation = triangles, mapper = mapper)
end

function Jutul.triangulate_mesh(mesh::EmbeddedReservoirMesh; outer = false, flatten = true, kwarg...)
    matrix = Jutul.triangulate_mesh(UnstructuredMesh(mesh.matrix_mesh); outer = outer, flatten = flatten, kwarg...)
    fractures = Jutul.triangulate_mesh(mesh.fracture_mesh; flatten = flatten, kwarg...)
    points = vcat(matrix.points, fractures.points)
    if flatten
        offset = size(matrix.points, 1)
        triangles = vcat(matrix.triangulation, fractures.triangulation .+ offset)
    else
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
