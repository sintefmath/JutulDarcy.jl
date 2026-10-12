"""
    UnstructuredMesh(mesh::EmbeddedReservoirMesh)

Convert a combined EFM or DFM mesh to a three-dimensional `UnstructuredMesh` for
plotting and geometric operations. Matrix cells and nodes retain their indices;
fracture cells follow the matrix cells in simulation order. Each fracture patch
is extruded symmetrically about its plane using its physical volume divided by
its area, including aperture overrides supplied during domain setup.

The original connection graph and face indices are retained. Added connections
are non-neighboring faces without nodes, while the extruded fracture surfaces
are boundary faces. The matrix geometry is copied without cutting it, so the
fracture prisms overlap the matrix. Use the unified domain's transmissibilities
and storage for simulation, rather than recomputing them from this geometric mesh.
"""
function Jutul.UnstructuredMesh(mesh::EmbeddedReservoirMesh)
    T = Jutul.float_type(mesh)
    matrix = UnstructuredMesh(mesh.matrix_mesh)
    fractures = mesh.fracture_mesh
    nm = number_of_cells(matrix)
    nc = number_of_cells(mesh)
    nf = number_of_faces(mesh)
    points = Point3{T}.(matrix.node_points)
    face_nodes = copy(matrix.faces.faces_to_nodes.vals)
    face_pos = copy(matrix.faces.faces_to_nodes.pos)
    # Added EFM/DFM connections are NNCs, rather than shared geometric faces.
    append!(face_pos, fill(last(face_pos), nf-number_of_faces(matrix)))
    N = tpfv_geometry(mesh).neighbors
    neighbors = [(N[1, f], N[2, f]) for f in 1:nf]
    cells_to_faces = embedded_cell_faces(N, nc)

    boundary_nodes = copy(matrix.boundary_faces.faces_to_nodes.vals)
    boundary_pos = copy(matrix.boundary_faces.faces_to_nodes.pos)
    boundary_cells = copy(matrix.boundary_faces.neighbors)
    cell_boundary_faces = copy(matrix.boundary_faces.cells_to_faces.vals)
    cell_boundary_pos = copy(matrix.boundary_faces.cells_to_faces.pos)
    parents = Int[]
    for i in 1:number_of_cells(fractures)
        vertices, faces = embedded_fracture_prism(mesh.matrix_mesh, fractures,
            tpfv_geometry(mesh).volumes[nm+i], i)
        offset = length(points)
        append!(points, vertices)
        for nodes in faces
            append!(boundary_nodes, nodes .+ offset)
            push!(boundary_pos, length(boundary_nodes)+1)
            push!(boundary_cells, nm+i)
            push!(cell_boundary_faces, length(boundary_cells))
        end
        push!(cell_boundary_pos, length(cell_boundary_faces)+1)
        push!(parents, fracture_parent_cell(fractures, mesh.matrix_mesh, i))
    end
    cell_map = collect(1:nm)
    if !isnothing(matrix.cell_map)
        cell_map = collect(matrix.cell_map)
    end
    append!(cell_map, cell_map[parents])
    converted = UnstructuredMesh(
        cells_to_faces,
        Jutul.IndirectionMap(cell_boundary_faces, cell_boundary_pos),
        Jutul.IndirectionMap(face_nodes, face_pos),
        Jutul.IndirectionMap(boundary_nodes, boundary_pos),
        points, neighbors, boundary_cells;
        structure = matrix.structure,
        cell_map = cell_map,
        z_is_depth = Jutul.mesh_z_is_depth(matrix))
    # Matrix fault and orientation tags still refer to the same entities.
    for (entity, tags) in pairs(Jutul.mesh_entity_tags(matrix).tags)
        for (group, values) in pairs(tags.tags)
            for (value, indices) in pairs(values)
                Jutul.set_mesh_entity_tag!(converted, entity, group, value, copy(indices))
            end
        end
    end
    return converted
end

function embedded_fracture_prism(matrix, fractures::EmbeddedFractureMesh, volume, cell)
    patch = fractures.geometry.patches[cell]
    frame = fractures.geometry.frames[patch.fracture]
    aperture = volume/patch.area
    return fracture_prism(patch.parts, frame, aperture, fractures.geometry.tolerance)
end

function embedded_fracture_prism(matrix, fractures::Jutul.EmbeddedMesh, volume, cell)
    T = eltype(matrix.node_points[1])
    if cell <= length(fractures.parent_faces)
        face = fractures.parent_faces[cell]
        points = matrix.node_points[matrix.faces.faces_to_nodes[face]]
        tol = T(64)*eps(T)*max(one(T), maximum(norm, points))
        frame = fracture_frame(PolygonFracture(points; float_type = T), tol)
        area, center = patch_geometry(frame.parts, frame)
        vertices, faces = fracture_prism(frame.parts, frame, volume/area, tol)
    else
        # Retained intersections are line cells. Give them a square cross section
        # whose area reproduces their physical (aperture squared) storage.
        mesh = UnstructuredMesh(fractures)
        face = first(mesh.faces.cells_to_faces[cell])
        nodes = mesh.faces.faces_to_nodes[face]
        a, b = Point3{T}.(mesh.node_points[nodes])
        segment_length = norm(b-a)
        aperture = sqrt(volume/segment_length)
        direction = (b-a)/segment_length
        u, normal = plane_basis(direction)
        frame = FractureFrame((a+b)/2, -normal, u, direction, Vector{Point2{T}}[])
        part = Point2{T}[(-aperture/2, -segment_length/2), (aperture/2, -segment_length/2),
            (aperture/2, segment_length/2), (-aperture/2, segment_length/2)]
        tol = T(64)*eps(T)*max(one(T), norm(frame.origin), segment_length)
        vertices, faces = fracture_prism([part], frame, aperture, tol)
    end
    return vertices, faces
end
