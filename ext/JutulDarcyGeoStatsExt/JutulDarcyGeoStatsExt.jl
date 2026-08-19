module JutulDarcyGeoStatsExt
    using JutulDarcy, Jutul, Random, LinearAlgebra, StaticArrays
    using GeoStats

    """
        JutulDarcy.generate_perm_poro(dims, box_lengths = (1.0, 1.0, 1.0); kwargs...)

    Generate porosity and permeability realizations over a box-shaped domain.

    The default `porosity_process` is a Gaussian process with spherical covariance,
    where the horizontal and vertical correlation lengths are scaled to the box
    dimensions. A custom generator can also be supplied for deterministic or
    problem-specific sampling.
    """
    function JutulDarcy.generate_perm_poro(
        dims::NTuple{3, Int},
        box_lengths::NTuple{3, <:Real} = (1.0, 1.0, 1.0);
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
        permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
        perm_from_poro = JutulDarcy.kozeny_carman_permeability,
    )
        nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
        all(>(0), dims) || throw(ArgumentError("All grid dimensions must be positive."))
        all(>(0), box_lengths) || throw(ArgumentError("All box lengths must be positive."))

        phimin, phimax = porosity_bounds
        0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))

        if !isnothing(seed)
            Random.seed!(seed)
        end

        box_origin_f = _to_float_tuple(box_origin)
        box_lengths_f = _to_float_tuple(box_lengths)

        spacing = ntuple(i -> Float64(box_lengths_f[i]) / dims[i], 3)
        volumes = fill(Float64(prod(spacing)), dims...)
        points = JutulDarcy.grid_points(dims, box_origin_f, spacing)

        if isnothing(porosity_process)
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
            porosity_process = GaussianProcess(
                SphericalCovariance(ranges = normalized_ranges),
                0.0,
            )
        end

        grid = CartesianGrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0); dims = dims)
        realizations = NamedTuple[]
        for _ in 1:nrealizations
            if porosity_process isa Function
                zporo = porosity_process(dims...)
            else
                zporo = _geostats_field_to_array(rand(porosity_process, grid).field, dims)
            end
            porosity = clamp.(porosity_mean .+ porosity_std .* zporo, phimin, phimax)
            try
                permeability = perm_from_poro(porosity; permeability_bounds = permeability_bounds)
            catch
                permeability = perm_from_poro(porosity; bounds = permeability_bounds)
            end
            push!(realizations, (
                porosity = porosity,
                permeability = permeability,
                points = points,
                volumes = volumes,
            ))
        end
        return nrealizations == 1 ? only(realizations) : realizations
    end

    """
        JutulDarcy.map_to_domain(mesh, tab::GeoStats.GeoTable; coordinate_mapping = nothing, kwargs...)

    Convert a GeoStats table sampled on a 3D grid into a reservoir property field by
    mapping the grid cell centroids to the enclosing reservoir cells.
    """
    function JutulDarcy.map_to_domain(mesh::JutulMesh, tab::GeoStats.GeoTable; coordinate_mapping = nothing, kwargs...)

        grid = tab.geometry
        paramdim(grid) == 3 || throw(ArgumentError("Expected a 3D GeoStats grid."))

        values = tab.field
        pts = centroid.(grid)
        pts = [SVector{3, Float64}([p.coords.x.val, p.coords.y.val, p.coords.z.val]) for p in pts]
        pts = coordinate_mapping === nothing ? pts : coordinate_mapping.(pts)

        property = (values = values, points = pts)

        return JutulDarcy.map_to_domain(mesh, property; kwargs...)

    end

end