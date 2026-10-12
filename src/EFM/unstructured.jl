"""
    UnstructuredMesh(mesh::EmbeddedReservoirMesh)

Convert a combined EFM mesh to a three-dimensional `UnstructuredMesh` for
plotting and geometric operations. Matrix cells and nodes retain their indices;
fracture cells follow the matrix cells in simulation order. Each fracture patch
is extruded symmetrically about its plane using its physical volume divided by
its area, including aperture overrides supplied during domain setup.

The original connection graph and face indices are retained. EFM connections
are non-neighboring faces without nodes, while the extruded fracture surfaces
are boundary faces. The matrix geometry is copied without cutting it, so the
fracture prisms overlap the matrix. Use the EFM domain's transmissibilities and
storage for simulation, rather than recomputing them from this geometric mesh.
"""
function Jutul.UnstructuredMesh(mesh::EmbeddedReservoirMesh)
    matrix = UnstructuredMesh(mesh.matrix_mesh)
    fractures = mesh.fracture_mesh
    nm = number_of_cells(matrix)
    nc = number_of_cells(mesh)
    nf = number_of_faces(mesh)
    points = Point3.(matrix.node_points)
    face_nodes = copy(matrix.faces.faces_to_nodes.vals)
    face_pos = copy(matrix.faces.faces_to_nodes.pos)
    # EFM connections are NNCs, rather than shared geometric faces.
    append!(face_pos, fill(last(face_pos), nf-number_of_faces(matrix)))
    N = mesh.geometry.neighbors
    neighbors = [(N[1, f], N[2, f]) for f in 1:nf]
    cells_to_faces = embedded_cell_faces(neighbors, nc)

    boundary_nodes = copy(matrix.boundary_faces.faces_to_nodes.vals)
    boundary_pos = copy(matrix.boundary_faces.faces_to_nodes.pos)
    boundary_cells = copy(matrix.boundary_faces.neighbors)
    cell_boundary_faces = copy(matrix.boundary_faces.cells_to_faces.vals)
    cell_boundary_pos = copy(matrix.boundary_faces.cells_to_faces.pos)
    parents = Int[]
    for (i, patch) in enumerate(fractures.patches)
        frame = fractures.frames[patch.fracture]
        aperture = mesh.geometry.volumes[nm+i]/patch.area
        vertices, faces = fracture_prism(patch.parts, frame, aperture, fractures.tolerance)
        offset = length(points)
        append!(points, vertices)
        for nodes in faces
            append!(boundary_nodes, nodes .+ offset)
            push!(boundary_pos, length(boundary_nodes)+1)
            push!(boundary_cells, nm+i)
            push!(cell_boundary_faces, length(boundary_cells))
        end
        push!(cell_boundary_pos, length(cell_boundary_faces)+1)
        push!(parents, first(patch.matrix_cells))
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

function embedded_cell_faces(neighbors, nc)
    counts = zeros(Int, nc+1)
    for (l, r) in neighbors
        counts[l+1] += 1
        counts[r+1] += 1
    end
    pos = cumsum(counts) .+ 1
    vals = Vector{Int}(undef, last(pos)-1)
    next = copy(pos)
    for (face, cells) in enumerate(neighbors)
        for cell in cells
            vals[next[cell]] = face
            next[cell] += 1
        end
    end
    return Jutul.IndirectionMap(vals, pos)
end

"""Extrude convex patch parts, removing walls inside a concave patch."""
function fracture_prism(parts, frame, aperture, tol)
    points = Point2[]
    polygons = [[point_index!(points, p, tol) for p in orient_polygon(part)] for part in parts]
    caps = Vector{Int}[]
    edges = Set{Tuple{Int, Int}}()
    for polygon in polygons
        cap = Int[]
        for i in eachindex(polygon)
            a, b = polygon[i], polygon[mod1(i+1, length(polygon))]
            origin, direction = points[a], points[b]-points[a]
            distance = norm(direction)
            direction /= distance
            # Split shared edges at all patch vertices. Clipping can produce
            # partial edge overlaps between the convex parts of one cell.
            splits = Tuple{Float64, Int}[(0.0, a), (distance, b)]
            for (j, p) in enumerate(points)
                t = dot(p-origin, direction)
                if tol < t < distance-tol && abs(cross2(direction, p-origin)) <= tol
                    push!(splits, (t, j))
                end
            end
            sort!(splits; by = first)
            for k in 1:length(splits)-1
                l, r = splits[k][2], splits[k+1][2]
                push!(cap, l)
                if (r, l) in edges
                    delete!(edges, (r, l))
                else
                    push!(edges, (l, r))
                end
            end
        end
        push!(caps, cap)
    end
    n = length(points)
    lower = [unproject(frame, p)-aperture/2*frame.normal for p in points]
    upper = [unproject(frame, p)+aperture/2*frame.normal for p in points]
    faces = Vector{Int}[]
    for cap in caps
        push!(faces, reverse(cap))
        push!(faces, cap .+ n)
    end
    for (l, r) in sort!(collect(edges))
        push!(faces, [l, r, r+n, l+n])
    end
    return vcat(lower, upper), faces
end
