module JutulDarcyGeoStatsExt
    using JutulDarcy, Jutul, Random, LinearAlgebra
    using GeoStats: CartesianGrid, GaussianProcess, SphericalCovariance

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
        box = (
            dims = dims,
            lengths = box_lengths_f,
            origin = box_origin_f,
        )

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
                box = box,
            ))
        end
        return realizations
    end

    """
        JutulDarcy.map_realization_to_reservoir_domain(domain, realization; copy_domain = true)

    Map a box-defined GeoStats realization onto an arbitrary 3D JutulDarcy
    reservoir domain.

    For each reservoir cell, all realization centroids inside the cell are
    aggregated with arithmetic porosity averaging and harmonic permeability
    averaging.
    """
    function JutulDarcy.map_realization_to_reservoir_domain(domain, realization; copy_domain::Bool = true)
        hasproperty(realization, :porosity) || throw(ArgumentError("realization must provide a porosity field."))
        hasproperty(realization, :permeability) || throw(ArgumentError("realization must provide a permeability field."))
        hasproperty(realization, :box) || throw(ArgumentError("realization must provide box metadata."))

        mesh = reservoir_mesh(domain)
        mesh_u = UnstructuredMesh(mesh)
        dim(mesh_u) == 3 || throw(ArgumentError("Expected a 3D reservoir mesh."))

        box = realization.box
        dims = box.dims
        porosity = realization.porosity
        permeability = realization.permeability
        size(porosity) == dims || throw(ArgumentError("Porosity field size $(size(porosity)) does not match box dimensions $dims."))
        size(permeability) == dims || throw(ArgumentError("Permeability field size $(size(permeability)) does not match box dimensions $dims."))

        axes_xyz = _box_centroid_axes(box)
        cache = _prepare_cell_mapping_cache(mesh_u)

        nc = number_of_cells(mesh_u)
        porosity_cells = Vector{Float64}(undef, nc)
        permeability_cells = Vector{Float64}(undef, nc)

        for cell in 1:nc
            bounds = cache.bounds[cell]
            ix = _axis_index_range(axes_xyz[1], bounds[1][1], bounds[2][1])
            iy = _axis_index_range(axes_xyz[2], bounds[1][2], bounds[2][2])
            iz = _axis_index_range(axes_xyz[3], bounds[1][3], bounds[2][3])

            sum_poro = 0.0
            sum_inv_perm = 0.0
            count = 0

            for k in iz, j in iy, i in ix
                pt = (axes_xyz[1][i], axes_xyz[2][j], axes_xyz[3][k])
                if _point_inside_cell(mesh_u, cache, cell, pt)
                    phi = porosity[i, j, k]
                    K = permeability[i, j, k]
                    K > 0 || throw(ArgumentError("Permeability must be strictly positive for harmonic averaging."))
                    sum_poro += phi
                    sum_inv_perm += inv(K)
                    count += 1
                end
            end

            count > 0 || throw(ArgumentError("Reservoir cell $cell did not capture any GeoStats cell centroids. Refine the realization grid or adjust the realization box."))
            porosity_cells[cell] = sum_poro / count
            permeability_cells[cell] = count / sum_inv_perm
        end

        mapped = copy_domain ? deepcopy(domain) : domain
        mapped[:porosity] = porosity_cells
        mapped[:permeability] = permeability_cells
        return mapped
    end

    function _to_float_tuple(values::NTuple{3, <:Real})
        return ntuple(i -> Float64(values[i]), 3)
    end

    function _geostats_field_to_array(field, dims::NTuple{3, Int})
        values = collect(field)
        length(values) == prod(dims) || throw(ArgumentError("GeoStats field has $(length(values)) values, expected $(prod(dims))."))
        return reshape(Float64.(vec(values)), dims)
    end

    function _box_centroid_axes(box)
        spacing = ntuple(i -> box.lengths[i] / box.dims[i], 3)
        return ntuple(i -> box.origin[i] .+ ((1:box.dims[i]) .- 0.5) .* spacing[i], 3)
    end

    function _prepare_cell_mapping_cache(mesh_u)
        geometry = tpfv_geometry(mesh_u)
        boundary_normals = zeros(size(geometry.boundary_centroids))
        for bface in axes(boundary_normals, 2)
            cell = mesh_u.boundary_faces.neighbors[bface]
            boundary_normals[:, bface] .= geometry.boundary_centroids[:, bface] .- geometry.cell_centroids[:, cell]
            boundary_normals[:, bface] ./= norm(boundary_normals[:, bface])
        end
        bounds = [JutulDarcy.cell_node_bounds(mesh_u, cell) for cell in 1:number_of_cells(mesh_u)]
        return (
            normals = geometry.normals,
            face_centroids = geometry.face_centroids,
            boundary_normals = boundary_normals,
            boundary_centroids = geometry.boundary_centroids,
            bounds = bounds,
        )
    end

    function _axis_index_range(axis, low, high)
        start = searchsortedfirst(axis, low)
        stop = searchsortedlast(axis, high)
        if start > length(axis) || stop < 1 || start > stop
            return 1:0
        end
        return start:stop
    end

    function _point_inside_cell(mesh_u, cache, cell, pt)
        for face in mesh_u.faces.cells_to_faces[cell]
            sgn = mesh_u.faces.neighbors[face][1] == cell ? 1.0 : -1.0
            if !_inside_half_space(pt, cache.normals, cache.face_centroids, face, sgn)
                return false
            end
        end
        for bface in mesh_u.boundary_faces.cells_to_faces[cell]
            if !_inside_half_space(pt, cache.boundary_normals, cache.boundary_centroids, bface, 1.0)
                return false
            end
        end
        return _point_in_bounds(pt, cache.bounds[cell])
    end

    function _inside_half_space(pt, normals, centroids, index, scale)
        val = 0.0
        for d in 1:3
            val += scale * normals[d, index] * (pt[d] - centroids[d, index])
        end
        return val <= 0.0
    end

    function _point_in_bounds(pt, bounds)
        low, high = bounds
        for d in 1:3
            if pt[d] < low[d] || pt[d] > high[d]
                return false
            end
        end
        return true
    end
end
