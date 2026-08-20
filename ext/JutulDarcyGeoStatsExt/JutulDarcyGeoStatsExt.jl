module JutulDarcyGeoStatsExt
    using JutulDarcy, Jutul, Random, LinearAlgebra, StaticArrays
    using GeoStats

    """
        generate_perm_poro(geometry::GeoStats.Meshes.GeometrySet; kwargs...)

    Generate correlated porosity and permeability realizations on a GeoStats
    geometry set.

    The geometry is first normalized to the unit cube so the Gaussian process is
    evaluated in a numerically stable coordinate system. The resulting fields are
    returned as a `Dict` with `:porosity` and `:permeability`, and optionally
    `:geo_table` when `extra_output = true`.
    """
    function generate_perm_poro(geometry::GeoStats.Meshes.GeometrySet;
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
        extra_output = false,)

        nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
        phimin, phimax = porosity_bounds
        0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))

        if !isnothing(seed)
            Random.seed!(seed)
        end

        transformed, points, scale = normalize_geometry_set(geometry)
        sx, sy, sz = scale

        if isnothing(porosity_process)
            horizontal_default = max(sx, sy) / 3
            vertical_default = sz / 3
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
            ranges = (
                horizontal_range / max(sx, eps(Float64)),
                horizontal_range / max(sy, eps(Float64)),
                vertical_range / max(sz, eps(Float64)),
            )
            porosity_process = GaussianProcess(
                SphericalCovariance(ranges = ranges),
                0.0,
            )
        end

        realizations = Dict[]
        for _ in 1:nrealizations
            tab = nothing
            if porosity_process isa Function
                zporo = porosity_process(length(points))
            else
                tab = rand(porosity_process, transformed)
                zporo = tab.field
            end
            porosity = clamp.(porosity_mean .+ porosity_std .* zporo, phimin, phimax)
            permeability = perm_from_poro(porosity; bounds = permeability_bounds)

            realization = Dict{Symbol, Any}(
                :porosity => porosity,
                :permeability => permeability,
            )

            if extra_output && !isnothing(tab)
                realization[:geo_table] = tab
            end
            push!(realizations, realization)
        end

        return nrealizations == 1 ? only(realizations) : realizations
    end

    """
        generate_perm_poro(mesh::JutulMesh; kwargs...)

    Generate a porosity/permeability realization on the cell-centroid geometry
    of a Jutul mesh.

    The mesh is converted to a `GeometrySet` of centroid points, normalized to
    the unit cube, sampled with the GeoStats-backed generator, and returned as a
    `Dict` with `:porosity` and `:permeability`.
    """
    function JutulDarcy.generate_perm_poro(mesh::JutulMesh; kwargs...)
        
        geometry = mesh_to_geometry_set(mesh)
        return generate_perm_poro(geometry; kwargs...)
        
    end

    """
        generate_perm_poro(domain::DataDomain; kwargs...)

    Generate a porosity/permeability realization on the cell-centroid geometry
    of a `DataDomain` and write the mapped values onto a copy of that domain.

    If `extra_output = true`, the returned tuple also includes the GeoStats table
    produced during sampling.
    """
    function JutulDarcy.generate_perm_poro(domain::DataDomain; kwargs...)
        
        realization = JutulDarcy.generate_perm_poro(physical_representation(domain); kwargs...)
        
        new_domain = deepcopy(domain)
        new_domain[:porosity] = realization[:porosity]
        new_domain[:permeability] = realization[:permeability]

        if haskey(realization, :geo_table)
            out = (new_domain, realization[:geo_table])
        else
            out = new_domain
        end

        return out
    end

    """
        mesh_to_geometry_set(mesh::JutulMesh)

    Convert a Jutul mesh to a GeoStats `GeometrySet` built from the cell
    centroids of the mesh.
    """
    function mesh_to_geometry_set(mesh::JutulMesh)

        geometry = tpfv_geometry(mesh)
        pts = [GeoStats.Meshes.Point(p[1], p[2], p[3]) for p in eachcol(geometry.cell_centroids)]

        return GeoStats.Meshes.GeometrySet(pts)

    end

    """
        normalize_geometry_set(geometry::GeoStats.Meshes.GeometrySet)

    Translate and scale the centroids of a `GeometrySet` to the unit cube.

    Returns the normalized geometry, the original centroid points, and the scale
    factors applied along each axis.
    """
    function normalize_geometry_set(geometry::GeoStats.Meshes.GeometrySet)
        points = centroid.(geometry)
        length(points) > 0 || throw(ArgumentError("GeometrySet must contain at least one point."))

        xs = [p.coords.x.val for p in points]
        ys = [p.coords.y.val for p in points]
        zs = [p.coords.z.val for p in points]
        xmin, xmax = minimum(xs), maximum(xs)
        ymin, ymax = minimum(ys), maximum(ys)
        zmin, zmax = minimum(zs), maximum(zs)

        sx = max(xmax - xmin, eps(Float64))
        sy = max(ymax - ymin, eps(Float64))
        sz = max(zmax - zmin, eps(Float64))

        transformed = GeoStats.Meshes.GeometrySet([
            GeoStats.Meshes.Point(
                (p.coords.x.val - xmin) / sx,
                (p.coords.y.val - ymin) / sy,
                (p.coords.z.val - zmin) / sz,
            ) for p in points
        ])

        return transformed, points, (sx, sy, sz)
    end

end
