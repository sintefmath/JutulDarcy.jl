"""
    discrete_fracture_mesh(matrix_mesh, fractures; min_cut_fraction=0.0, kwarg...)

Cut a matrix mesh using the same fracture objects as `embedded_fracture_mesh`.
Returns `(matrix_mesh, fracture_mesh, info)`, where the fracture mesh is a DFM
`Jutul.EmbeddedMesh`. `info[:cell_index]` maps cut cells to original cells;
`:face_index` and `:boundary_face_index` provide the corresponding face maps,
and `:fracture_id` maps DFM fracture cells to the input objects.

Cuts along polygon edges make finite fracture boundaries conforming, including
concave polygons and the polygonal approximation of an ellipse. Auxiliary cuts
are restricted to original cells intersecting the fracture. They refine the
matrix without creating additional fractures. The input mesh is unchanged.
"""
function discrete_fracture_mesh(matrix_mesh::JutulMesh, fractures;
        tolerance = nothing,
        min_cut_fraction = 0.0,
        merge_faces = false,
        intersection_strategy = :star_delta
    )
    intersection_strategy in (:star_delta, :remove) || throw(ArgumentError("intersection_strategy must be :star_delta or :remove"))
    !merge_faces || throw(ArgumentError("merge_faces must be false to preserve cut-face mappings"))
    original = UnstructuredMesh(matrix_mesh)
    any(f -> isempty(original.faces.faces_to_nodes[f]), 1:number_of_faces(original)) &&
        throw(ArgumentError("DFM cutting requires a matrix mesh without non-neighboring faces"))
    embedded = embedded_fracture_mesh(matrix_mesh, fractures; tolerance = tolerance)
    geometry = embedded.geometry
    tol = geometry.tolerance
    cut = original
    cell_index = collect(1:number_of_cells(original))
    face_index = collect(1:number_of_faces(original))
    boundary_index = collect(1:number_of_boundary_faces(original))
    function apply_cut(plane, parents, frame)
        cells = findall(c -> cell_index[c] in parents &&
            fracture_intersects_cell(cut, c, frame, tol), eachindex(cell_index))
        cut, step = Jutul.cut_mesh(cut, plane; cells = cells, extra_out = true,
            min_cut_fraction = min_cut_fraction, merge_faces = merge_faces)
        cell_index = cell_index[step[:cell_index]]
        face_index = compose_fracture_mapping(face_index, step[:face_index])
        boundary_index = compose_fracture_mapping(boundary_index, step[:boundary_face_index])
        return nothing
    end
    for (id, object) in enumerate(embedded.fractures)
        parents = Set{Int}()
        for patch in geometry.patches
            if patch.fracture == id
                union!(parents, patch.matrix_cells)
            end
        end
        frame = geometry.frames[id]
        points = fracture_polygon(object)
        for i in eachindex(points)
            a, b = points[i], points[mod1(i+1, length(points))]
            normal = normalize(cross(b-a, frame.normal))
            apply_cut(Jutul.PlaneCut(a, normal), parents, frame)
        end
        apply_cut(Jutul.PlaneCut(frame.origin, frame.normal), parents, frame)
    end
    ids = zeros(Int, number_of_faces(cut))
    for (id, frame) in enumerate(geometry.frames)
        for face in eachindex(ids)
            nodes = cut.faces.faces_to_nodes[face]
            points = cut.node_points[nodes]
            all(p -> abs(dot(p-frame.origin, frame.normal)) <= 8tol, points) || continue
            parts = convex_parts([project(frame, p) for p in points], tol)
            area, _ = patch_geometry(parts, frame)
            covered = 0.0
            for part in parts, boundary in frame.parts
                clipped = clip_polygon(part, boundary, tol)
                if length(clipped) >= 3
                    covered += polygon_area_centroid(clipped)[1]
                end
            end
            covered > tol^2 || continue
            allowance = max(16tol*sqrt(area), 128eps(typeof(tol))*area)
            area-covered <= allowance || throw(ArgumentError("DFM face $face is only partially covered by fracture $id; reduce min_cut_fraction to zero"))
            ids[face] == 0 || throw(ArgumentError("coplanar overlapping fractures must be merged before cutting"))
            ids[face] = id
        end
    end
    faces = findall(>(0), ids)
    fracture_ids = ids[faces]
    for id in eachindex(embedded.fractures)
        id in fracture_ids || throw(ArgumentError("fracture $id produced no interior DFM faces; boundary-only fractures are not supported"))
    end
    cell_map = cell_index
    if !isnothing(original.cell_map)
        cell_map = original.cell_map[cell_index]
    end
    cut = UnstructuredMesh(cut.faces.cells_to_faces, cut.boundary_faces.cells_to_faces,
        cut.faces.faces_to_nodes, cut.boundary_faces.faces_to_nodes,
        cut.node_points, cut.faces.neighbors, cut.boundary_faces.neighbors;
        structure = original.structure, cell_map = cell_map,
        z_is_depth = Jutul.mesh_z_is_depth(original))
    for (entity, mapping) in ((Cells(), cell_index), (Faces(), face_index),
            (BoundaryFaces(), boundary_index))
        tags = Jutul.mesh_entity_tags(original)[entity]
        for (group, values) in pairs(tags.tags), (value, indices) in pairs(values)
            selected = Set(indices)
            mapped = findall(parent -> parent in selected, mapping)
            Jutul.set_mesh_entity_tag!(cut, entity, group, value, mapped)
        end
        if entity != Cells() && haskey(tags, :orientation)
            vertical, horizontal = Int[], Int[]
            for face in findall(==(0), mapping)
                normal = Jutul.face_normal(cut, face, entity)
                if abs(normal[3]) < max(abs(normal[1]), abs(normal[2]))
                    push!(vertical, face)
                else
                    push!(horizontal, face)
                end
            end
            Jutul.set_mesh_entity_tag!(cut, entity, :orientation, :vertical, vertical)
            Jutul.set_mesh_entity_tag!(cut, entity, :orientation, :horizontal, horizontal)
        end
    end
    fracture_mesh = Jutul.EmbeddedMesh(cut, faces; intersection_strategy = intersection_strategy)
    info = Dict{Symbol, Any}(:cell_index => cell_index, :face_index => face_index,
        :boundary_face_index => boundary_index, :fracture_id => fracture_ids,
        :fracture_faces => faces)
    return cut, fracture_mesh, info
end

function fracture_intersects_cell(mesh, cell, frame, tol)
    sections, _ = cell_section(mesh, cell, frame.origin, frame.normal, frame.u, frame.v, tol)
    for section in sections, part in convex_parts(section, tol; resolve_winding = true), boundary in frame.parts
        clipped = clip_polygon(part, boundary, tol)
        if length(clipped) >= 3 && polygon_area_centroid(clipped)[1] > tol^2
            return true
        end
    end
    return false
end

function compose_fracture_mapping(previous, current)
    out = zeros(Int, length(current))
    for (i, parent) in enumerate(current)
        if parent > 0
            out[i] = previous[parent]
        end
    end
    return out
end

"""
    discrete_fracture_domain(matrix, fractures; aperture=nothing, permeability=nothing, porosity=nothing, kwarg...)

Create cut-cell matrix and DFM fracture domains from EFM fracture objects.
Returns `(matrix=..., fractures=..., info=...)`. Matrix cell properties are
inherited through `info[:cell_index]`; fracture properties use the input object
IDs and accept the same per-object overrides as `embedded_fracture_domain`.
Geometry keywords are forwarded to `discrete_fracture_mesh`.
"""
function discrete_fracture_domain(matrix::DataDomain, objects;
        aperture = nothing,
        hydraulic_aperture = nothing,
        permeability = nothing,
        porosity = nothing,
        kwarg...
    )
    objects = FractureCollection(objects)
    cut, mesh, info = discrete_fracture_mesh(physical_representation(matrix), objects; kwarg...)
    cut_matrix = cut_fracture_matrix_domain(matrix, cut, info)
    ids = info[:fracture_id]
    defaults = fracture_property_defaults(objects, Jutul.float_type(cut))
    a = fracture_property(aperture, defaults.aperture[ids], ids, "aperture")
    k = fracture_property(permeability, defaults.permeability[ids], ids, "permeability")
    phi = fracture_property(porosity, defaults.porosity[ids], ids, "porosity")
    ah = a
    if !isnothing(hydraulic_aperture)
        ah = fracture_property(hydraulic_aperture, a, ids, "hydraulic_aperture")
        all(>=(0), ah) || throw(ArgumentError("hydraulic_aperture must be non-negative"))
    end
    if isnothing(permeability) && (!isnothing(aperture) || !isnothing(hydraulic_aperture))
        k = ah.^2/12
    end
    all(>(0), a) || throw(ArgumentError("aperture must be positive"))
    all(>=(0), k) || throw(ArgumentError("permeability must be non-negative"))
    all(x -> 0 < x <= 1, phi) || throw(ArgumentError("porosity must lie in (0, 1]"))
    fractures = JutulDarcy.fracture_domain(mesh, cut_matrix; aperture = a, permeability = k, porosity = phi)
    # Project tensor matrix permeability onto each fracture's normal, just as
    # in EFM exchange, rather than indexing a tensor as a scalar array.
    fmc = FractureMatrixConnection()
    mc, fc = fractures[:matrix_cells, fmc], fractures[:connection_cells, fmc]
    normals = fractures[:cell_normals]
    fractures[:matrix_permeability, fmc] = [normal_conductivity(cut_matrix[:permeability], cell,
        normals[:, fc[i]]) for (i, cell) in enumerate(mc)]
    fractures[:fracture_id, Cells()] = ids
    return (matrix = cut_matrix, fractures = fractures, info = info)
end

function cut_fracture_matrix_domain(matrix, mesh, info)
    domain = reservoir_domain(mesh)
    geometry_keys = keys(DataDomain(mesh).data)
    for (key, (value, entity)) in pairs(matrix)
        key in geometry_keys && continue
        if entity == Cells()
            domain[key, entity] = fracture_cell_subset(value, info[:cell_index])
        elseif entity == Faces() || entity == BoundaryFaces()
            mapping = info[:face_index]
            if entity == BoundaryFaces()
                mapping = info[:boundary_face_index]
            end
            if key == :transmissibility_multiplier
                mapped = ones(eltype(value), length(mapping))
            elseif key == :transmissibility_override
                mapped = fill(eltype(value)(NaN), length(mapping))
            elseif key == :transmissibilities
                # Explicit coefficients are transferred after cell properties
                # and multipliers, so newly introduced faces can use TPFA.
                continue
            else
                throw(ArgumentError("cannot automatically transfer face property $key to a cut mesh; apply it after discrete_fracture_domain"))
            end
            for (i, parent) in enumerate(mapping)
                if parent > 0
                    mapped[i] = value[parent]
                    if key == :transmissibility_override && isfinite(mapped[i])
                        mapped[i] *= domain[:areas][i]/matrix[:areas][parent]
                    end
                end
            end
            domain[key, entity] = mapped
        elseif isnothing(entity) || entity == Jutul.NoEntity()
            domain[key, entity] = deepcopy(value)
        end
    end
    if haskey(matrix, :transmissibilities, Faces())
        values = JutulDarcy.reservoir_transmissibility(domain)
        for (face, parent) in enumerate(info[:face_index])
            if parent > 0
                values[face] = matrix[:transmissibilities][parent]*domain[:areas][face]/matrix[:areas][parent]
            end
        end
        domain[:transmissibilities, Faces()] = values
    end
    return domain
end

fracture_cell_subset(values::AbstractVector, cells) = values[cells]
fracture_cell_subset(values::AbstractMatrix, cells) = values[:, cells]
