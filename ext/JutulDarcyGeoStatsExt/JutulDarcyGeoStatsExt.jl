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
        box_lengths::NTuple{3, <:Real} = (1.0, 1.0, 1.0),
        box_origin::NTuple{3, <:Real} = (0.0, 0.0, 0.0);
        nrealizations::Int = 1,
        seed = nothing,
        porosity_process = nothing,
        porosity_mean::Real = 0.20,
        porosity_std::Real = 0.05,
        porosity_bounds::Tuple{<:Real, <:Real} = (0.05, 0.95),
        correlation_range = nothing,
        horizontal_correlation_range = nothing,
        vertical_correlation_range = nothing,
        permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
        perm_from_poro = JutulDarcy.kozeny_carman_permeability,
        extra_output = false,
    )
        nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
        all(>(0), dims) || throw(ArgumentError("All grid dimensions must be positive."))
        all(>(0), box_lengths) || throw(ArgumentError("All box lengths must be positive."))

        phimin, phimax = porosity_bounds
        0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))

        if !isnothing(seed)
            Random.seed!(seed)
        end

        to_float_tuple(t::NTuple{3, <:Real}) = (Float64(t[1]), Float64(t[2]), Float64(t[3]))

        box_origin_f = to_float_tuple(box_origin)
        box_lengths_f = to_float_tuple(box_lengths)

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
        # Map from unit cube to box domain
        points = centroid.(grid)
        points = [box_origin_f .+ SVector{3, Float64}([
            p.coords.x.val * box_lengths_f[1],
            p.coords.y.val * box_lengths_f[2],
            p.coords.z.val * box_lengths_f[3]
            ]) for p in points]

        # Scale volumes to physical box dimensions
        volumes = measure.(grid)
        volumes = [v.val * prod(box_lengths_f) for v in volumes]

        realizations = Dict[]
        for _ in 1:nrealizations
            if porosity_process isa Function
                zporo = porosity_process(dims...)
            else
                tab = rand(porosity_process, grid)
                zporo = tab.field
            end
            porosity = clamp.(porosity_mean .+ porosity_std .* zporo, phimin, phimax)
            permeability = perm_from_poro(porosity; bounds = permeability_bounds)
            
            realization = Dict{Symbol, Any}(
                :porosity => (values = porosity, points = points, volumes = volumes),
                :permeability => (values = permeability, points = points, volumes = volumes)
            )

            if extra_output
                realization[:geo_table] = tab
            end
            push!(realizations, realization)

        end

        return nrealizations == 1 ? only(realizations) : realizations

    end

end