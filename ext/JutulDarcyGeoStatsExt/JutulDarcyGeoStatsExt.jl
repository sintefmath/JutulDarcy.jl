module JutulDarcyGeoStatsExt
    using JutulDarcy, Jutul, Random, LinearAlgebra, StaticArrays
    using GeoStats

    """
        JutulDarcy.kozeny_carman_permeability(porosity; kozeny_constant = 1.0*si_unit(:darcy),
            permeability_bounds = (1e-20, Inf))

    Compute scalar permeability from porosity using a Kozeny-Carman style relation.
    """
    function JutulDarcy.kozeny_carman_permeability(
        porosity;
        kozeny_constant::Real = 1.0*si_unit(:darcy),
        permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
    )
        kozeny_constant > 0 || throw(ArgumentError("kozeny_constant must be positive."))

        kmin, kmax = permeability_bounds
        0.0 < kmin <= kmax || throw(ArgumentError("permeability_bounds must satisfy 0 < min <= max."))

        numerator = kozeny_constant .* porosity.^3
        denominator = (1 .- porosity).^2
        permeability = numerator ./ denominator
        return clamp.(permeability, kmin, kmax)
    end

    """
        JutulDarcy.setup_perm_poro_realizations(dims, box_lengths; kwargs...)

    Generate GeoStats-based porosity realizations over a box domain and derive
    correlated permeability realizations from porosity.

    The GeoStats process is sampled on a unit cube internally. Correlation ranges
    are specified in physical units and normalized by box lengths before creating
    the covariance model.

    The returned realization stores physical box metadata in `realization.box`,
    consumed by [`map_realization_to_reservoir_domain`](@ref).
    """
    function JutulDarcy.setup_perm_poro_realizations(
        dims::NTuple{3, Int},
        box_lengths::NTuple{3, <:Real};
        nrealizations::Int = 1,
        seed = nothing,
        box_origin::NTuple{3, <:Real} = (0.0, 0.0, 0.0),
        porosity_process = nothing,
        porosity_mean::Real = 0.20,
        porosity_std::Real = 0.05,
        porosity_bounds::Tuple{<:Real, <:Real} = (0.05, 0.95),
        correlation_range = nothing,
        horizontal_correlation_range = nothing,
        vertical_correlation_range = nothing,
        grain_diameter::Real = 1e-5,
        permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
    )
        nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
        all(>(0), dims) || throw(ArgumentError("All grid dimensions must be positive."))
        all(>(0), box_lengths) || throw(ArgumentError("All box lengths must be positive."))

        _ = grain_diameter # reserved for future constitutive variants

        phimin, phimax = porosity_bounds
        0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))

        if !isnothing(seed)
            Random.seed!(seed)
        end

        box_origin_f = _to_float_tuple(box_origin)
        box_lengths_f = _to_float_tuple(box_lengths)

        horizontal_default = max(box_lengths_f[1], box_lengths_f[2]) / 3
        vertical_default = box_lengths_f[3] / 3
        if isnothing(correlation_range)
            horizontal_range = isnothing(horizontal_correlation_range) ? horizontal_default : Float64(horizontal_correlation_range)
            vertical_range = isnothing(vertical_correlation_range) ? vertical_default : Float64(vertical_correlation_range)
        else
            base_range = Float64(correlation_range)
            horizontal_range = isnothing(horizontal_correlation_range) ? base_range : Float64(horizontal_correlation_range)
            vertical_range = isnothing(vertical_correlation_range) ? base_range : Float64(vertical_correlation_range)
        end
        horizontal_range > 0 || throw(ArgumentError("horizontal_correlation_range must be positive."))
        vertical_range > 0 || throw(ArgumentError("vertical_correlation_range must be positive."))

        normalized_ranges = (
            horizontal_range / box_lengths_f[1],
            horizontal_range / box_lengths_f[2],
            vertical_range / box_lengths_f[3],
        )

        grid = CartesianGrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0); dims = dims)
        if isnothing(porosity_process)
            porosity_process = GaussianProcess(
                SphericalCovariance(ranges = normalized_ranges),
                0.0,
            )
        end

        porosity_latent = rand(porosity_process, grid, nrealizations)

        realizations = NamedTuple[]
        for i in 1:nrealizations
            zporo = _geostats_field_to_array(porosity_latent[i].field, dims)
            porosity = clamp.(porosity_mean .+ porosity_std .* zporo, phimin, phimax)
            permeability = JutulDarcy.kozeny_carman_permeability(
                porosity;
                permeability_bounds = permeability_bounds,
            )
            push!(realizations, (
                porosity = porosity,
                permeability = permeability,
                grid = grid,
            ))
        end
        return realizations
    end

    function JutulDarcy.map_to_domain!(domain, property, grid, name; kwargs...)
        values = map_to_domain(domain, property; name = name, kwargs...)
        domain[name] = values
        return domain
    end

    function JutulDarcy.map_to_domain!(domain, properties::Dict, grid; kwargs...)
        for (name, property) in properties
            values = map_to_domain(domain, property; name = name, kwargs...)
            domain[name] = values
        end
        return domain
    end

    function JutulDarcy.map_to_domain(mesh::JutulMesh, tab::GeoStats.GeoTable, coordinate_mapping = nothing; kwargs...)

        grid = tab.geometry
        paramdim(grid) == 3 || throw(ArgumentError("Expected a 3D GeoStats grid."))

        values = tab.field
        pts = centroid.(grid)
        pts = [SVector{3, Float64}([p.coords.x.val, p.coords.y.val, p.coords.z.val]) for p in pts]
        pts = coordinate_mapping === nothing ? pts : coordinate_mapping.(pts)

        return JutulDarcy.map_to_domain(mesh, values, pts; kwargs...)

    end

    function JutulDarcy.map_to_domain(mesh::JutulMesh, values, points; mapping::Symbol = :mean, info_level = 0, name = missing)

        if points isa AbstractMatrix
            nrows, ncols = size(points)
            nrows in (2, 3) || throw(ArgumentError("Point matrix must have 2 or 3 rows, got $(nrows)."))
            pts = [SVector{nrows, Float64}(ntuple(i -> Float64(points[i, j]), nrows)) for j in 1:ncols]
        else
            pts = [SVector{length(pt), Float64}(Float64.(pt)) for pt in points]
        end

        total_input_cells = length(values)
        total_input_cells == length(pts) || throw(ArgumentError("Property length $(length(values)) does not match the number of input points $(length(pts))."))

        mesh_u = UnstructuredMesh(mesh)
        dim(mesh_u) == 3 || throw(ArgumentError("Expected a 3D reservoir mesh."))

        geometry = tpfv_geometry(mesh_u)
        cell_lookup = create_cell_lookup(mesh_u, geometry, pts)
        mapped_values, total_assigned = aggregate_property_to_cells(values, cell_lookup, mapping)

        ignored_input_cells = total_input_cells - total_assigned
        if info_level > 0
            approx_cells_per_output = total_input_cells / max(number_of_cells(mesh_u), 1)
            str = ifelse(ismissing(name), "", "$name: ")
            @info str*"Mapping statistics: approx. $(approx_cells_per_output) input cells per output cell; $(ignored_input_cells) input cells were outside the mesh and ignored."
        end

        return mapped_values
        
    end

    function JutulDarcy.map_to_domain(domain::DataDomain, property; kwargs...)

        mesh = domain.mesh
        return map_to_domain(mesh, property; kwargs...)

    end

    function aggregate_property_to_cells(property_array, cell_lookup, mapping::Symbol)
        nc = maximum(cell_lookup)
        mapped_values = fill(NaN, nc)
        counts = zeros(Int, nc)
        sum_values = zeros(Float64, nc)
        sum_inv_values = zeros(Float64, nc)
        total_assigned = 0

        for idx in eachindex(property_array)
            cell = cell_lookup[idx]
            cell == 0 && continue
            total_assigned += 1
            value = Float64(property_array[idx])
            counts[cell] += 1
            if mapping === :mean
                sum_values[cell] += value
            elseif mapping === :harmonic_mean
                value > 0 || throw(ArgumentError("Property values must be strictly positive for harmonic averaging."))
                sum_inv_values[cell] += inv(value)
            else
                throw(ArgumentError("mapping must be :mean or :harmonic_mean."))
            end
        end

        for cell in 1:nc
            count = counts[cell]
            count == 0 && continue
            if mapping === :mean
                mapped_values[cell] = sum_values[cell] / count
            else
                mapped_values[cell] = count / sum_inv_values[cell]
            end
        end

        return mapped_values, total_assigned
    end

    function create_cell_lookup(mesh, geometry, pts)
        total = length(pts)
        cell_lookup = zeros(Int, total)

        for idx in 1:total
            pt = pts[idx]
            cells = Jutul.find_enclosing_cells(
                mesh,
                [pt, pt];
                geometry = geometry,
                n = 1,
                limit_box = false,
                cells = 1:number_of_cells(mesh),
            )
            if isempty(cells)
                cell_lookup[idx] = 0
            else
                cell_lookup[idx] = only(cells)
            end
        end
        return cell_lookup
    end

end

    # function _to_float_tuple(values::NTuple{3, <:Real})
    #     return ntuple(i -> Float64(values[i]), 3)
    # end