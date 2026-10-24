using Test, Jutul, JutulDarcy, GeoEnergyIO, GeoStats

function dipping_cpg_mesh(; nx = 6, ny = 5, nz = 3)
    xs = collect(range(0.0, 100.0, length = nx + 1))
    ys = collect(range(0.0, 100.0, length = ny + 1))
    top = [10.0 + 0.1x + 0.05y for x in xs, y in ys]
    horizons = [top .+ 4.0*k for k in 0:nz]
    grid = cpgrid_from_horizons(xs, ys, horizons)
    return mesh_from_grid_section(grid)
end

@testset "map_property_to_cells" begin
    source = CartesianMesh((20, 20, 8), (100.0, 100.0, 20.0))
    sgeo = tpfv_geometry(source)
    x_source = sgeo.cell_centroids[1, :]
    @testset "upscaling and downscaling" begin
        coarse = CartesianMesh((5, 5, 2), (100.0, 100.0, 20.0))
        cgeo = tpfv_geometry(coarse)
        # Aligned upscaling of a linear field is exact
        up = map_property_to_cells(coarse, source, x_source)
        @test up ≈ cgeo.cell_centroids[1, :]
        # Downscaling gives each fine cell the value of the coarse cell containing it
        x_coarse = cgeo.cell_centroids[1, :]
        down = map_property_to_cells(source, coarse, x_coarse)
        @test down ≈ 10.0 .+ 20.0 .* floor.(x_source ./ 20.0)
        # Matrix values are mapped row by row
        K = vcat(x_source', 2 .* x_source', 3 .* x_source')
        Kc = map_property_to_cells(coarse, source, K)
        @test size(Kc) == (3, number_of_cells(coarse))
        @test Kc[2, :] ≈ 2 .* up
    end
    @testset "non-uniform target meshes" begin
        dx = 10.0 .^ range(-0.5, 1.2, length = 12)
        dx = 100.0 .* dx ./ sum(dx)
        tensor = CartesianMesh((12, 5, 2), (dx, fill(20.0, 5), fill(10.0, 2)))
        cpg = dipping_cpg_mesh()
        for mesh in (tensor, cpg)
            geo = tpfv_geometry(mesh)
            vals, W = map_property_to_cells(mesh, source, x_source, extra_out = true, outside = :nearest)
            # Overlap volumes add up to the cell volumes
            @test vec(sum(W, dims = 2)) ≈ geo.volumes rtol = 1e-8
            @test all(isfinite, vals)
            # Constant fields are preserved for all averages
            for mapping in (:arithmetic, :geometric, :harmonic)
                c = map_property_to_cells(mesh, source, fill(0.25, length(x_source)), mapping = mapping, weights = W)
                @test all(isapprox.(c, 0.25))
            end
        end
        # Volume-weighted mapping preserves the total for a mesh covering the
        # source exactly, also when the mesh is finer than the source in places
        vals = map_property_to_cells(tensor, source, x_source)
        @test sum(vals .* tpfv_geometry(tensor).volumes) ≈ sum(x_source .* sgeo.volumes) rtol = 1e-3
        # Ordering of averages
        perm = exp.(x_source ./ 20.0)
        a = map_property_to_cells(tensor, source, perm, mapping = :arithmetic)
        g = map_property_to_cells(tensor, source, perm, mapping = :geometric)
        h = map_property_to_cells(tensor, source, perm, mapping = :harmonic)
        @test all(a .>= g .- 1e-12)
        @test all(g .>= h .- 1e-12)
    end
    @testset "mesh outside source" begin
        big = CartesianMesh((4, 4, 2), (200.0, 100.0, 20.0))
        v = map_property_to_cells(big, source, x_source, outside = :nearest)
        @test all(isfinite, v)
        @test v[4] ≈ 97.5
        v = map_property_to_cells(big, source, x_source, outside = :ignore)
        @test all(isnan, v[[3, 4]])
        @test all(isfinite, v[[1, 2]])
        @test_throws ArgumentError map_property_to_cells(big, source, x_source, outside = :error)
        @test_throws ArgumentError map_property_to_cells(big, source, x_source[1:10])
        @test_throws ArgumentError map_property_to_cells(big, source, x_source, mapping = :mean)
    end
    @testset "scattered points" begin
        coarse = CartesianMesh((5, 5, 2), (100.0, 100.0, 20.0))
        cgeo = tpfv_geometry(coarse)
        pts = sgeo.cell_centroids
        v, lookup = map_property_to_cells(coarse, pts, x_source, volumes = sgeo.volumes, extra_out = true)
        @test v ≈ cgeo.cell_centroids[1, :]
        @test all(lookup .> 0)
        # Coarse points to a fine mesh use the nearest point
        fine_v = map_property_to_cells(source, cgeo.cell_centroids, cgeo.cell_centroids[1, :])
        @test fine_v ≈ 10.0 .+ 20.0 .* floor.(x_source ./ 20.0)
        @test_throws ArgumentError map_property_to_cells(coarse, pts[1:2, :], x_source)
        # Points outside the mesh are ignored
        small = CartesianMesh((2, 2, 2), (50.0, 50.0, 10.0))
        v, lookup = map_property_to_cells(small, pts, x_source, extra_out = true)
        @test count(iszero, lookup) == length(x_source) - 10*10*4
        @test v ≈ [12.5, 37.5, 12.5, 37.5, 12.5, 37.5, 12.5, 37.5]
    end
    @testset "map_property_to_cells!" begin
        coarse = CartesianMesh((5, 5, 2), (100.0, 100.0, 20.0))
        domain = reservoir_domain(coarse)
        perm = exp.(x_source ./ 20.0)
        props = Dict(:porosity => fill(0.3, length(x_source)), :permeability => perm)
        map_property_to_cells!(domain, source, props)
        @test all(domain[:porosity] .≈ 0.3)
        @test domain[:permeability] ≈ map_property_to_cells(coarse, source, perm, mapping = :geometric)
        map_property_to_cells!(domain, source, props, mapping = :harmonic)
        @test domain[:permeability] ≈ map_property_to_cells(coarse, source, perm, mapping = :harmonic)
        map_property_to_cells!(domain, sgeo.cell_centroids, perm, :permeability)
        @test domain[:permeability] ≈ map_property_to_cells(coarse, sgeo.cell_centroids, perm)
        # GeoStats grids can be used as source
        grid = CartesianGrid((0.0, 0.0, 0.0), (100.0, 100.0, 20.0), dims = (20, 20, 8))
        @test map_property_to_cells(coarse, grid, x_source) ≈ map_property_to_cells(coarse, source, x_source)
    end
end

@testset "generate_perm_poro" begin
    g = CartesianMesh((10, 10, 4), (100.0, 100.0, 20.0))
    nc = number_of_cells(g)
    r = generate_perm_poro(g, seed = 1)
    @test r isa AbstractDict
    @test length(r[:porosity]) == nc
    @test length(r[:permeability]) == nc
    @test all(0.05 .<= r[:porosity] .<= 0.95)
    @test r[:permeability] ≈ kozeny_carman_permeability(r[:porosity])
    @test r[:porosity] == generate_perm_poro(g, seed = 1)[:porosity]
    @test r[:porosity] == generate_perm_poro(g, rng = Xoshiro(1))[:porosity]
    @test r[:porosity] != generate_perm_poro(g, seed = 2)[:porosity]
    # The IJK grid of a Cartesian mesh is the mesh itself
    grid = CartesianGrid((0.0, 0.0, 0.0), (100.0, 100.0, 20.0), dims = (10, 10, 4))
    @test r[:porosity] ≈ generate_perm_poro(grid, seed = 1)[:porosity]

    rs = generate_perm_poro(g, nrealizations = 3, seed = 3)
    @test rs isa AbstractVector
    @test length(rs) == 3
    @test rs[1][:porosity] != rs[2][:porosity]

    for method in (:ijk, :centroids, :box)
        rm = generate_perm_poro(g, method = method, seed = 4, extra_out = true)
        @test length(rm[:porosity]) == nc
        @test all(0.05 .<= rm[:porosity] .<= 0.95)
        @test haskey(rm, :geotable)
    end
    @test_throws ArgumentError generate_perm_poro(g, method = :nearest)

    # Anisotropic permeability
    rk = generate_perm_poro(g, seed = 5, kv_kh = 0.1)
    @test size(rk[:permeability]) == (3, nc)
    @test rk[:permeability][3, :] ≈ 0.1 .* rk[:permeability][1, :]

    # Custom processes and transforms
    rf = generate_perm_poro(g, porosity_process = (rng, n) -> randn(rng, n), seed = 6)
    @test length(rf[:porosity]) == nc
    rt = generate_perm_poro(g, seed = 7, porosity_transform = z -> 0.3, porosity_bounds = (0.0, 0.4))
    @test all(rt[:porosity] .== 0.3)

    # Corner-point mesh with dipping layers
    cpg = dipping_cpg_mesh()
    for method in (:ijk, :centroids, :box)
        rc = generate_perm_poro(cpg, method = method, seed = 8)
        @test length(rc[:porosity]) == number_of_cells(cpg)
    end

    # Mesh with inactive cells
    actnum = ones(Int, 10, 20, 2)
    actnum[1:5, 1:10, 1] .= 0
    m = reservoir_mesh((10, 20, 2), (100.0, 200.0, 10.0), actnum = actnum)
    @test length(generate_perm_poro(m, seed = 9)[:porosity]) == number_of_cells(m)

    # DataDomain
    domain = reservoir_domain(g)
    rd = generate_perm_poro(domain, nrealizations = 2, seed = 10)
    @test length(rd) == 2
    rd = generate_perm_poro!(domain, seed = 11)
    @test domain[:porosity] == rd[:porosity]
    @test domain[:permeability] == rd[:permeability]
    @test_throws ArgumentError generate_perm_poro!(domain, nrealizations = 2)
end
