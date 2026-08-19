module JutulDarcyGeoStatsExt
    using JutulDarcy, Jutul, Random, LinearAlgebra, StaticArrays
    using GeoStats

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
    function JutulDarcy.generate_perm_poro_realizations(
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