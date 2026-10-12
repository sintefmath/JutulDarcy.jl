function well_segment_fracture_hits(mesh::EmbeddedFractureMesh{T}, index, a, b) where T
    tol = mesh.geometry.tolerance
    delta = b-a
    segment_length = norm(delta)
    hits = Tuple{T, Int, Point3{T}}[]
    if segment_length > tol
        for cell in query_cells(index, min.(a, b), max.(a, b), tol)
            patch = mesh.geometry.patches[cell]
            frame = mesh.geometry.frames[patch.fracture]
            da, db = dot(a-frame.origin, frame.normal), dot(b-frame.origin, frame.normal)
            if abs(da-db) > tol
                fraction = da/(da-db)
                if -tol/segment_length <= fraction <= 1+tol/segment_length
                    fraction = clamp(fraction, zero(T), one(T))
                    point = a+fraction*delta
                    projected = project(frame, point)
                    inside = any(part -> all(i -> cross2(part[mod1(i+1, length(part))]-part[i],
                        projected-part[i]) >= -tol*norm(part[mod1(i+1, length(part))]-part[i]),
                        eachindex(part)), patch.parts)
                    if inside
                        push!(hits, (fraction, cell, point))
                    end
                end
            elseif abs(da) <= tol && abs(db) <= tol
                # A well segment in the fracture plane gets one completion for
                # each positive-length portion inside a clipped fracture cell.
                pa, pb = project(frame, a), project(frame, b)
                for (lo, hi) in line_intervals(patch.parts, pa, (pb-pa)/segment_length,
                        tol; limits = (zero(T), segment_length))
                    fraction = (lo+hi)/(2segment_length)
                    push!(hits, (fraction, cell, a+fraction*delta))
                end
            end
        end
        sort!(hits; by = hit -> (hit[1], hit[2]))
        # Shared patch edges must not create two completions for one fracture.
        seen = Tuple{Int, T}[]
        filter!(hits) do hit
            fracture = mesh.geometry.patches[hit[2]].fracture
            keep = !any(pair -> pair[1] == fracture && abs(pair[2]-hit[1])*segment_length <= tol, seen)
            if keep
                push!(seen, (fracture, hit[1]))
            end
            keep
        end
    end
    return hits
end

function extend_embedded_well(well, fractures, index)
    mesh = physical_representation(fractures)
    w = physical_representation(well)
    T = Jutul.float_type(mesh)
    centers = well[:cell_centroids, Cells()]
    neighbors = NTuple{2, Int}[]
    old_indices, fracture_cells, well_cells = Int[], Int[], Int[]
    points = Point3{T}[]
    seen = Set{Tuple{Int, Int}}()
    nc = w.num_nodes
    for (segment, (left, right)) in enumerate(eachcol(w.neighborship))
        a, b = Point3{T}(centers[:, left]), Point3{T}(centers[:, right])
        hits = well_segment_fracture_hits(mesh, index, a, b)
        previous = left
        last_fraction = -one(T)
        for (fraction, cell, point) in hits
            if norm(point-a) <= mesh.geometry.tolerance
                node = left
            elseif norm(point-b) <= mesh.geometry.tolerance
                node = right
            elseif abs(fraction-last_fraction)*norm(b-a) <= mesh.geometry.tolerance
                node = previous
            else
                nc += 1
                node = nc
                push!(neighbors, (previous, node))
                push!(old_indices, segment)
            end
            if !((cell, node) in seen)
                push!(fracture_cells, cell)
                push!(well_cells, node)
                push!(points, point)
                push!(seen, (cell, node))
            end
            previous = node
            last_fraction = fraction
        end
        if previous != right
            push!(neighbors, (previous, right))
            push!(old_indices, segment)
        end
    end
    extended = well
    if !isempty(fracture_cells)
        N = reshape(reinterpret(Int, neighbors), 2, :)
        extended = JutulDarcy.extend_fractured_well(copy_fracture_data(well), fractures,
            N, fracture_cells, well_cells, old_indices)
        for (node, point) in zip(well_cells, points)
            extended[:cell_centroids, Cells()][:, node] .= point
        end
    end
    return extended
end

function extend_embedded_wells(wells, fractures)
    geometry = physical_representation(fractures).geometry
    index = fracture_patch_spatial_index(geometry.patches, geometry.frames)
    extended = DataDomain[]
    for well in wells
        w = physical_representation(well)
        if w isa JutulDarcy.MultiSegmentWell
            push!(extended, extend_embedded_well(well, fractures, index))
        else
            @warn "Well $(w.name) is a SimpleWell, which is currently not supported for extension with fracture perforations. The well will not be directly connected to any fractures."
            push!(extended, well)
        end
    end
    return extended
end

"""
Gravity parameter on `FracturePerforations`, shared by coupled EFM and DFM.
For each completion, returns `g*(z_fracture - z_well_node)` in m²/s² using that
completion's `self_fracture` node. The type has no fields. Unified models use
the ordinary perforation gravity parameter after merging completion lists.
"""
struct FracturePerforationGravityDifference <: ScalarVariable end
Jutul.associated_entity(::FracturePerforationGravityDifference) = JutulDarcy.FracturePerforations()

function Jutul.default_parameter_values(domain, model, ::FracturePerforationGravityDifference, name)
    well = physical_representation(domain)
    cells = well.perforations.self_fracture
    centers = domain[:cell_centroids, Cells()]
    perforations = domain[:perforation_centroids_frac, JutulDarcy.FracturePerforations()]
    differences = JutulDarcy.gravity_constant.*(perforations[3, :] .- centers[3, cells])
    return differences
end
