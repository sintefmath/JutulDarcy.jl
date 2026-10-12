@testset "Independent EFM fracture grids" begin
    matrix = reservoir_domain(reservoir_mesh((4, 4, 2), (4.0, 4.0, 2.0)))
    original = deepcopy(matrix.data)
    base = RectangleFracture([2, 2, 1], [0, 0, 1], (4, 4); aperture = 0.01)
    reference = embedded_fracture_domain(matrix, base)
    function matrix_exchange_totals(domain)
        geometry = physical_representation(domain).geometry
        totals = zeros(number_of_cells(matrix))
        for (cell, area) in zip(geometry.matrix_cells, geometry.matrix_areas)
            totals[cell] += area
        end
        return totals
    end
    reference_volumes = combined_fracture_domain(matrix, reference; subtract_fracture_volume = true)[:volumes][1:32]
    for resolution in ((1, 1), (2, 2), (8, 8))
        fracture = RectangleFracture([2, 2, 1], [0, 0, 1], (4, 4);
            aperture = 0.01, resolution = resolution)
        domain = embedded_fracture_domain(matrix, fracture)
        mesh = physical_representation(domain)
        grid = only(mesh.geometry.grids)
        nx, ny = resolution
        @test number_of_cells(domain) == nx*ny
        @test number_of_faces(domain) == (nx-1)*ny+nx*(ny-1)
        @test grid.resolution == resolution
        @test grid.spacing == (4/nx, 4/ny)
        @test sum(domain[:volumes]) ≈ sum(reference[:volumes])
        @test matrix_exchange_totals(domain) ≈ matrix_exchange_totals(reference)
        @test domain[:matrix_face_coverage] ≈ reference[:matrix_face_coverage]
        @test all(p -> sum(p.matrix_weights) ≈ 1, mesh.geometry.patches)
        @test all(p -> length(p.matrix_cells) == length(p.matrix_weights) == length(p.matrix_areas), mesh.geometry.patches)
        unified = combined_fracture_domain(matrix, domain; subtract_fracture_volume = true)
        @test unified[:volumes][1:32] ≈ reference_volumes
        converted = UnstructuredMesh(physical_representation(combined_fracture_domain(matrix, domain)))
        for cell in 33:number_of_cells(converted)
            center, volume = Jutul.compute_centroid_and_measure(converted, Cells(), cell)
            @test center ≈ domain[:cell_centroids][:, cell-32]
            @test volume ≈ domain[:volumes][cell-32]
        end
        @test isequal(matrix.data, original)
    end

    # Unequal matrix overlaps must retain their own storage/exchange fractions.
    uneven = reservoir_domain(UnstructuredMesh(CartesianMesh((2, 1, 1), ([1.0, 2.0], 1.0, 1.0))))
    f = RectangleFracture([1.5, 0.5, 0.5], [0, 0, 1], (3, 1); aperture = 0.01, resolution = 1)
    domain = embedded_fracture_domain(uneven, f)
    patch = only(physical_representation(domain).geometry.patches)
    @test patch.matrix_cells == [1, 2]
    @test patch.matrix_weights ≈ [1/3, 2/3]
    @test patch.matrix_areas ≈ [2, 4]
    @test combined_fracture_domain(uneven, domain; subtract_fracture_volume = true)[:volumes][1:2] ≈ [0.99, 1.98]

    # Subdivide inside one matrix cell and check the analytical planar TPFA stencil.
    unit_matrix = reservoir_domain(reservoir_mesh((1, 1, 1), (4.0, 4.0, 2.0)))
    rectangle = RectangleFracture([2, 2, 1], [0, 0, 1], (2, 3);
        resolution = (2, 3), angle = 0.3, aperture = 0.01, permeability = 1e-12)
    domain = embedded_fracture_domain(unit_matrix, rectangle)
    mesh = physical_representation(domain)
    @test number_of_cells(domain) == 6
    @test all(h -> h ≈ 1, only(mesh.geometry.grids).spacing)
    for (face, (l, r)) in enumerate(eachcol(domain[:neighbors]))
        distance = norm(domain[:cell_centroids][:, l]-domain[:cell_centroids][:, r])
        @test domain[:transmissibilities][face] ≈ 1e-14*mesh.geometry.finite_volume.areas[face]/distance
    end

    for T in (Float32, Float64, BigFloat)
        fracture = RectangleFracture(T.([2, 2, 1]), [0, 0, 1], (4, 4); cell_length = (1.5, 1.0))
        d = embedded_fracture_domain(matrix, fracture; float_type = T)
        grid = only(physical_representation(d).geometry.grids)
        @test grid isa EFM.FractureGrid{T}
        @test grid.resolution == (3, 4)
        @test grid.spacing[1] <= T(1.5) && grid.spacing[2] <= one(T)
        @test sum(d[:volumes]) ≈ T(16)*fracture.properties.aperture rtol = 1e-5
    end
end

@testset "Per-fracture controls, polygons and disconnected regions" begin
    matrix = reservoir_domain(reservoir_mesh((4, 4, 2), (4.0, 4.0, 2.0)))
    objects = Any[
        RectangleFracture([2, 2, 0.5], [0, 0, 1], (4, 4); resolution = (2, 3)),
        EllipticalFracture([2, 2, 1.5], [0, 0, 1], (1.5, 1.0); cell_length = 0.6)]
    domain = embedded_fracture_domain(matrix, objects;
        resolution = [nothing, (5, 4)], cell_length = [0.5, nothing])
    geometry = physical_representation(domain).geometry
    @test geometry.grids[1].resolution == (8, 8)
    @test geometry.grids[2].resolution == (5, 4)
    reset = embedded_fracture_domain(matrix, objects; resolution = :matrix)
    @test all(isnothing, physical_representation(reset).geometry.grids)
    for id in 1:2
        reference = embedded_fracture_domain(matrix, objects[id])
        @test sum(domain[:volumes][domain[:fracture_id] .== id]) ≈ sum(reference[:volumes])
    end
    polygon = PolygonFracture([[0, 0, 0.5], [2, 0, 0.5], [2, 1, 0.5],
        [1, 1, 0.5], [1, 2, 0.5], [0, 2, 0.5]]; cell_length = (0.25, 0.4))
    concave = embedded_fracture_domain(matrix, polygon)
    @test sum(concave[:volumes])/polygon.properties.aperture ≈ 3
    @test all(isfinite, concave[:transmissibilities])
    inactive = reservoir_domain(extract_submesh(reservoir_mesh((3, 1, 1), (3.0, 1.0, 1.0)), [1, 3]))
    disconnected = RectangleFracture([1.5, 0.5, 0.5], [0, 0, 1], (3, 1); resolution = 1)
    d = embedded_fracture_domain(inactive, disconnected)
    @test number_of_cells(d) == 2
    @test number_of_faces(d) == 0
    @test all(p -> length(p.matrix_cells) == 1, physical_representation(d).geometry.patches)
    dfm_default = discrete_fracture_domain(matrix, RectangleFracture([2, 2, 1], [0, 0, 1], (4, 4)))
    dfm_configured = discrete_fracture_domain(matrix, RectangleFracture([2, 2, 1], [0, 0, 1], (4, 4); resolution = (8, 8)))
    @test dfm_default.info[:cell_index] == dfm_configured.info[:cell_index]
    @test dfm_default.info[:fracture_faces] == dfm_configured.info[:fracture_faces]
    @test dfm_default.fractures[:volumes] == dfm_configured.fractures[:volumes]

    for options in ((resolution = 0,), (resolution = (2, 0),), (resolution = (2, 1.5),),
            (cell_length = -1,), (cell_length = Inf,), (resolution = 2, cell_length = 1))
        @test_throws ArgumentError RectangleFracture([0, 0, 0], [0, 0, 1], (1, 1); options...)
        @test_throws ArgumentError embedded_fracture_domain(matrix, objects; options...)
    end
    @test_throws ArgumentError embedded_fracture_domain(matrix, objects; resolution = [2])
    @test_throws ArgumentError embedded_fracture_domain(matrix, objects; cell_length = [1])
    @test_throws ArgumentError embedded_fracture_domain(matrix, objects; resolution = 10^10)
end

@testset "Internally meshed EFM: coupled and unified flow" begin
    for system in (ImmiscibleSystem((AqueousPhase(), LiquidPhase())), :geothermal)
        matrix = reservoir_domain(reservoir_mesh((3, 3, 3), (3.0, 3.0, 3.0)); permeability = 1e-2si_unit(:darcy))
        fractures = [RectangleFracture([1.5, 1.5, 1.5], [1, 0, 0], (3, 3);
                resolution = (2, 2), aperture = 0.005, permeability = si_unit(:darcy)),
            RectangleFracture([1.5, 1.5, 1.5], [0, 1, 0], (3, 3);
                cell_length = 1.0, aperture = 0.005, permeability = si_unit(:darcy))]
        wells = [setup_well(matrix, [1]; name = :Injector, simple_well = true),
            setup_well(matrix, [27]; name = :Producer, simple_well = true)]
        models = [setup_fractured_reservoir_model(matrix, fractures, system;
            method = :efm, unified = unified, wells = wells) for unified in (false, true)]
        initial = (Pressure = 50si_unit(:bar),)
        injection = (density = 1000.0,)
        mix = [1.0]
        key = :Temperature
        if system == :geothermal
            initial = (initial..., Temperature = 363.15)
            injection = (injection..., temperature = 283.15)
        else
            initial = (initial..., Saturations = [0.0, 1.0])
            mix = [1.0, 0.0]
            key = :Saturations
        end
        control = Dict(:Injector => InjectorControl(TotalRateTarget(1e-5), mix; injection...),
            :Producer => ProducerControl(BottomHolePressureTarget(10si_unit(:bar))))
        results = []
        for model in models
            state0 = setup_reservoir_state(model; initial...)
            forces = setup_reservoir_forces(model; control = control)
            case = JutulCase(model, fill(20.0, 4), forces; state0 = state0)
            solver = GenericKrylov(preconditioner = ILUZeroPreconditioner(), relative_tolerance = 1e-10)
            result = simulate_reservoir(case; info_level = -1, max_nonlinear_iterations = 40,
                linear_solver = solver, tol_cnv = 1e-10, tol_cnve = 1e-10, timesteps = :none,
                tol_cnv_well = 1e-6, tol_cnve_well = 1e-6, max_timestep_cuts = 0)
            @test all(r -> length(r[:ministeps]) == 1 && only(r[:ministeps])[:success], result.result.reports)
            push!(results, result)
        end
        for (coupled, unified) in zip(results[1].result.states, results[2].result.states)
            @test unified[:Reservoir][:Pressure] ≈ vcat(coupled[:Reservoir][:Pressure], coupled[:Fractures][:Pressure]) rtol = 1e-7
            if key == :Saturations
                expected = hcat(coupled[:Reservoir][key], coupled[:Fractures][key])
            else
                expected = vcat(coupled[:Reservoir][key], coupled[:Fractures][key])
            end
            @test unified[:Reservoir][key] ≈ expected rtol = 1e-7 atol = 1e-7
        end
    end
end
