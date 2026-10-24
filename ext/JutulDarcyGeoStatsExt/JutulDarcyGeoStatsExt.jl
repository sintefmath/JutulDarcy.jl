module JutulDarcyGeoStatsExt
    using JutulDarcy, Jutul, StaticArrays
    using GeoStats
    import GeoStats: Random

    function JutulDarcy.generate_perm_poro(geometry::GeoStats.Domain;
        nrealizations::Int = 1,
        rng = Random.default_rng(),
        seed = nothing,
        porosity_process = nothing,
        porosity_mean::Real = 0.20,
        porosity_std::Real = 0.05,
        porosity_transform = nothing,
        porosity_bounds::Tuple{<:Real, <:Real} = (0.05, 0.95),
        correlation_range = nothing,
        horizontal_correlation_range = nothing,
        vertical_correlation_range = nothing,
        simulation_method = nothing,
        permeability_bounds::Tuple{<:Real, <:Real} = (1e-20, Inf),
        perm_from_poro = JutulDarcy.kozeny_carman_permeability,
        kv_kh = nothing,
        info_level = 0,
        extra_out = false
        )

        nrealizations > 0 || throw(ArgumentError("nrealizations must be positive."))
        phimin, phimax = porosity_bounds
        0.0 <= phimin < phimax < 1.0 || throw(ArgumentError("porosity_bounds must satisfy 0 <= min < max < 1."))
        isnothing(kv_kh) || kv_kh > 0 || throw(ArgumentError("kv_kh must be positive."))
        if !isnothing(seed)
            rng = Random.Xoshiro(seed)
        end
        if isnothing(porosity_process)
            porosity_process = default_porosity_process(geometry,
                correlation_range,
                horizontal_correlation_range,
                vertical_correlation_range
            )
        end
        if isnothing(porosity_transform)
            porosity_transform = z -> porosity_mean + porosity_std*z
        end

        n = nelements(geometry)
        if porosity_process isa Function
            tabs = fill(nothing, nrealizations)
        else
            # Simulate all realizations at once so that preprocessing (e.g.
            # factorization of the covariance matrix) is only done once.
            ensemble = rand(rng, porosity_process, geometry, nrealizations,
                method = simulation_method,
                showprogress = info_level > 0
            )
            tabs = [ensemble[i] for i in 1:nrealizations]
        end
        realizations = Dict{Symbol, Any}[]
        for tab in tabs
            if isnothing(tab)
                z = porosity_process(rng, n)
            else
                z = getproperty(tab, only(setdiff(propertynames(tab), (:geometry,))))
            end
            length(z) == n || throw(ArgumentError("porosity_process returned $(length(z)) values, expected $n."))
            porosity = porosity_transform.(z)
            nclamped = count(x -> x < phimin || x > phimax, porosity)
            porosity = clamp.(porosity, phimin, phimax)
            if info_level > 0
                jutul_message("GeoStats", "$(round(100*nclamped/n, sigdigits = 3))% of the porosity values were clamped to $porosity_bounds.")
            end
            permeability = perm_from_poro(porosity; bounds = permeability_bounds)
            if !isnothing(kv_kh)
                permeability = vcat(permeability', permeability', kv_kh.*permeability')
            end
            realization = Dict{Symbol, Any}(
                :porosity => porosity,
                :permeability => permeability,
            )
            if extra_out && !isnothing(tab)
                realization[:geotable] = tab
            end
            push!(realizations, realization)
        end
        return nrealizations == 1 ? only(realizations) : realizations
    end

    function JutulDarcy.generate_perm_poro(mesh::JutulMesh;
            method = JutulDarcy.mesh_has_ijk(mesh) ? :ijk : :centroids,
            ijk_spacing = missing,
            box_dims = missing,
            perm_mapping = :geometric,
            info_level = 0,
            kwarg...
        )
        dim(mesh) == 3 || throw(ArgumentError("Expected a 3D mesh."))
        if method == :ijk
            JutulDarcy.mesh_has_ijk(mesh) || throw(ArgumentError("method = :ijk requires a mesh with IJK structure (e.g. Cartesian or corner-point)."))
            dims = grid_dims_ijk(mesh)
            if ismissing(ijk_spacing)
                ijk_spacing = JutulDarcy.average_cell_spacing(mesh)
            end
            low, = JutulDarcy.mesh_bounding_box(mesh)
            grid = GeoStats.CartesianGrid(Tuple(low), Tuple(low) .+ dims.*ijk_spacing, dims = dims)
            out = JutulDarcy.generate_perm_poro(grid; info_level = info_level, kwarg...)
            ix = JutulDarcy.cell_ijk_linear_indices(mesh)
            for realization in realizations_vector(out)
                realization[:porosity] = realization[:porosity][ix]
                realization[:permeability] = select_cells(realization[:permeability], ix)
            end
        elseif method == :centroids
            geometry = tpfv_geometry(mesh)
            pts = [GeoStats.Point(p[1], p[2], p[3]) for p in eachcol(geometry.cell_centroids)]
            out = JutulDarcy.generate_perm_poro(GeoStats.GeometrySet(pts); info_level = info_level, kwarg...)
        elseif method == :box
            low, high = JutulDarcy.mesh_bounding_box(mesh)
            if ismissing(box_dims)
                box_dims = default_box_dims(mesh, low, high)
            end
            grid = GeoStats.CartesianGrid(Tuple(low), Tuple(high), dims = Tuple(box_dims))
            out = JutulDarcy.generate_perm_poro(grid; info_level = info_level, kwarg...)
            source = cartesian_grid_to_mesh(grid)
            weights = missing
            for realization in realizations_vector(out)
                realization[:porosity], weights = map_property_to_cells(mesh, source, realization[:porosity],
                    mapping = :arithmetic,
                    weights = weights,
                    info_level = info_level,
                    name = :porosity,
                    extra_out = true
                )
                realization[:permeability] = map_property_to_cells(mesh, source, realization[:permeability],
                    mapping = perm_mapping,
                    weights = weights
                )
            end
        else
            throw(ArgumentError("method must be :ijk, :centroids or :box, was :$method."))
        end
        return out
    end

    function JutulDarcy.generate_perm_poro(domain::DataDomain; kwarg...)
        return JutulDarcy.generate_perm_poro(physical_representation(domain); kwarg...)
    end

    function JutulDarcy.generate_perm_poro!(domain::DataDomain; nrealizations::Int = 1, kwarg...)
        nrealizations == 1 || throw(ArgumentError("generate_perm_poro! generates a single realization, use generate_perm_poro for more."))
        realization = JutulDarcy.generate_perm_poro(domain; kwarg...)
        domain[:porosity] = realization[:porosity]
        domain[:permeability] = realization[:permeability]
        return realization
    end

    """
        map_property_to_cells(mesh, grid::CartesianGrid, values; kwarg...)

    Map values defined on the elements of a GeoStats `CartesianGrid` (e.g. a
    realization from `rand(process, grid)`) to the cells of `mesh`. See the
    method taking a `CartesianMesh` source for details and keyword arguments.
    """
    function JutulDarcy.map_property_to_cells(mesh::JutulMesh, grid::GeoStats.CartesianGrid, values; kwarg...)
        return map_property_to_cells(mesh, cartesian_grid_to_mesh(grid), values; kwarg...)
    end

    function cartesian_grid_to_mesh(grid::GeoStats.CartesianGrid)
        length(size(grid)) == 3 || throw(ArgumentError("Expected a 3D CartesianGrid."))
        origin = ustrip.(GeoStats.to(minimum(grid)))
        spacing = ustrip.(GeoStats.spacing(grid))
        dims = size(grid)
        return CartesianMesh(dims, dims.*spacing, origin = collect(origin))
    end

    function default_porosity_process(geometry, correlation_range, horizontal_correlation_range, vertical_correlation_range)
        sx, sy, sz = ustrip.(GeoStats.sides(GeoStats.boundingbox(geometry)))
        if isnothing(correlation_range)
            horizontal_default = max(sx, sy)/3
            vertical_default = sz/3
        else
            horizontal_default = vertical_default = Float64(correlation_range)
        end
        horizontal_range = Float64(something(horizontal_correlation_range, horizontal_default))
        vertical_range = Float64(something(vertical_correlation_range, vertical_default))
        # Degenerate extents (e.g. a single layer of centroids) give zero range
        horizontal_range = max(horizontal_range, eps())
        vertical_range = max(vertical_range, eps())
        horizontal_range > 0 || throw(ArgumentError("horizontal_correlation_range must be positive."))
        vertical_range > 0 || throw(ArgumentError("vertical_correlation_range must be positive."))
        ranges = (horizontal_range, horizontal_range, vertical_range)
        return GaussianProcess(SphericalCovariance(ranges = ranges), 0.0)
    end

    function default_box_dims(mesh, low, high)
        # Resolve the average cell size of the mesh in each direction
        spacing = JutulDarcy.average_cell_spacing(mesh)
        extent = high - low
        return Tuple(max.(1, round.(Int, extent./spacing)))
    end

    realizations_vector(r::AbstractDict) = (r,)
    realizations_vector(r::AbstractVector) = r

    select_cells(k::AbstractVector, ix) = k[ix]
    select_cells(k::AbstractMatrix, ix) = k[:, ix]
end
