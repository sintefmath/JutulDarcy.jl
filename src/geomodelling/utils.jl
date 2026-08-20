"""
    map_to_domain(mesh, property::NamedTuple; mapping = :mean, cell_lookup = missing,
        info_level = 0, name = missing, extra_out = false)

Map a property defined on a set of input points to the cells of a Jutul reservoir
mesh.

The input property is expected to contain at least `:points` and `:values`.
Optionally, `:volumes` may be supplied for volume-weighted averaging. Supported
mappings are:
- `:mean` for arithmetic averaging
- `:harmonic_mean` for harmonic averaging
- `:volume_weighted_mean` for weighted arithmetic averaging using `volumes`

If `cell_lookup` is omitted, enclosing cells are computed from the mesh geometry.
The `info_level` keyword emits a summary log with the approximate input/output cell
ratio and the number of discarded points outside the mesh.
"""
function map_to_domain(mesh::JutulMesh, property::NamedTuple;
    mapping::Symbol = :mean,
    cell_lookup = missing,
    info_level = 0,
    name = missing,
    extra_out = false
    )

    points = property.points
    values = property.values
    volumes = get(property, :volumes, nothing)

    if points isa AbstractMatrix
        nrows, ncols = size(points)
        nrows in (2, 3) || throw(ArgumentError("Point matrix must have 2 or 3 rows, got $(nrows)."))
        pts = [SVector{nrows, Float64}(ntuple(i -> Float64(points[i, j]), nrows)) for j in 1:ncols]
    else
        pts = [SVector{length(pt), Float64}(Float64.(pt)) for pt in points]
    end

    total_input_cells = length(values)
    total_input_cells == length(pts) || throw(ArgumentError("Property length $(length(values)) does not match the number of input points $(length(pts))."))

    if !isnothing(volumes)
        total_input_cells == length(volumes) || throw(ArgumentError("Property length $(length(values)) does not match the number of volumes $(length(volumes))."))
    end

    mesh_u = UnstructuredMesh(mesh)
    dim(mesh_u) == 3 || throw(ArgumentError("Expected a 3D reservoir mesh."))

    if ismissing(cell_lookup)
        geometry = tpfv_geometry(mesh_u)
        cell_lookup = create_cell_lookup(mesh_u, geometry, pts)
    end
    mapped_values, total_assigned = aggregate_property_to_cells(values, cell_lookup, mapping; volumes = volumes)

    ignored_input_cells = total_input_cells - total_assigned
    if info_level > 0
        approx_cells_per_output = total_input_cells / max(number_of_cells(mesh_u), 1)
        str = ifelse(ismissing(name), "", "$name: ")
        @info str*"Mapping statistics: approx. $(approx_cells_per_output) \
        input cells per output cell; $(ignored_input_cells) input cells \
        were outside the mesh and ignored."
    end

    if extra_out
        return mapped_values, cell_lookup, ignored_input_cells
    else
        return mapped_values
    end
    
end

"""
    map_to_domain(domain, property::NamedTuple; kwargs...)

Convenience wrapper that maps a property defined on an input point set onto the
reservoir mesh represented by a `DataDomain`.
"""
function JutulDarcy.map_to_domain(domain::DataDomain, property::NamedTuple; kwargs...)

    mesh = physical_representation(domain)
    return map_to_domain(mesh, property; kwargs...)

end

"""
    map_to_domain!(domain, property::NamedTuple, name; kwargs...)

Mutating variant of `map_to_domain` that stores the mapped values on the domain
under the given property name.
"""
function map_to_domain!(domain::DataDomain, property::NamedTuple, name; kwargs...)

    values = map_to_domain(domain, property; name = name, kwargs...)
    domain[name] = values
    return domain

end

function map_to_domain!(domain::DataDomain, properties::Dict{Symbol, Any};
    cell_lookup = missing, kwargs...)

    for (name, property) in properties
        property isa NamedTuple || continue
        values, cell_lookup, _ = map_to_domain(domain, property;
            name = name, cell_lookup = cell_lookup, extra_out = true, kwargs...)
        domain[name] = values
    end

    return domain

end

function aggregate_property_to_cells(property_array, cell_lookup, mapping::Symbol; volumes = nothing)
    nc = maximum(cell_lookup)
    mapped_values = fill(NaN, nc)
    counts = zeros(Int, nc)
    sum_values = zeros(Float64, nc)
    sum_inv_values = zeros(Float64, nc)
    weighted_sum = zeros(Float64, nc)
    total_volume = zeros(Float64, nc)
    total_assigned = 0

    for idx in eachindex(property_array)
        cell = cell_lookup[idx]
        cell == 0 && continue
        total_assigned += 1
        value = Float64(property_array[idx])
        counts[cell] += 1
        if mapping === :mean || mapping === :volume_weighted_mean
            sum_values[cell] += value
        elseif mapping === :harmonic_mean
            value > 0 || throw(ArgumentError("Property values must be strictly positive for harmonic averaging."))
            sum_inv_values[cell] += inv(value)
        else
            throw(ArgumentError("mapping must be :mean, :volume_weighted_mean or :harmonic_mean."))
        end

        if !isnothing(volumes)
            vol = Float64(volumes[idx])
            vol > 0 || throw(ArgumentError("Volumes must be strictly positive for volume-weighted aggregation."))
            if mapping === :volume_weighted_mean
                weighted_sum[cell] += value * vol
                total_volume[cell] += vol
            end
        end
    end

    for cell in 1:nc
        count = counts[cell]
        count == 0 && continue
        if mapping === :mean
            mapped_values[cell] = sum_values[cell] / count
        elseif mapping === :harmonic_mean
            mapped_values[cell] = count / sum_inv_values[cell]
        elseif mapping === :volume_weighted_mean
            if total_volume[cell] > 0
                mapped_values[cell] = weighted_sum[cell] / total_volume[cell]
            else
                mapped_values[cell] = sum_values[cell] / count
            end
        end
    end

    return mapped_values, total_assigned
end

function create_cell_lookup(mesh, geometry, pts)
    total = length(pts)
    cell_lookup = zeros(Int, total)

    T = eltype(pts)
    normals = vec(reinterpret(T, geometry.normals))
    face_centroids = vec(reinterpret(T, geometry.face_centroids))
    boundary_centroids = vec(reinterpret(T, geometry.boundary_centroids))
    boundary_normals = vec(reinterpret(T, geometry.boundary_normals))
    
    covered = falses(number_of_cells(mesh))
    for idx in 1:total
        pt = pts[idx]
        cell = Jutul.find_enclosing_cell(mesh, pt, normals, face_centroids, boundary_normals, boundary_centroids)
        if isnothing(cell)
            cell_lookup[idx] = 0
        else
            cell_lookup[idx] = cell
            covered[cell] = true
        end
    end

    @assert all(covered) "Not all cells were covered by the input points. Consider increasing the number of input points or adjusting their distribution."

    return cell_lookup
end
