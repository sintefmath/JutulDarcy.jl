@testset "EFM well intersection geometry" begin
    matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0)))
    fracture = RectangleFracture([1, 1, 1], [1, 0, 0], (2, 2))
    domain = embedded_fracture_domain(matrix, fracture)
    centers = [0.5 1.5; 0.5 0.5; 0.2 1.6]
    well = setup_well(matrix, [1, 2]; simple_well = false, well_cell_centers = centers)
    before = deepcopy(well.data)
    system = SinglePhaseSystem(AqueousPhase())
    coupled, pc = setup_embedded_fracture_model(matrix, domain, system;
        wells = [well], extend_wells = true, extra_out = true)
    w = physical_representation(coupled[:Well].data_domain)
    @test length(w.perforations.fracture) == 1
    @test w.num_nodes == 3
    @test w.num_segments == 2
    point = coupled[:Well].data_domain[:cell_centroids][:, only(w.perforations.self_fracture)]
    @test point ≈ [1.0, 0.5, 0.9]
    fc = only(w.perforations.fracture)
    expected_gravity = JutulDarcy.gravity_constant*(domain[:cell_centroids][3, fc]-point[3])
    @test only(pc[:Well][:FracturePerforationGravityDifference]) ≈ expected_gravity
    @test abs(expected_gravity) > 1
    @test isequal(well.data, before)
    unified, pu = setup_fractured_reservoir_model(matrix, Any[fracture], system;
        method = :efm, unified = true, wells = [well], extend_wells = true, extra_out = true)
    uw = physical_representation(unified[:Well].data_domain)
    @test uw.perforations.reservoir == vcat(w.perforations.reservoir, w.perforations.fracture .+ number_of_cells(matrix))
    @test last(pu[:Well][:PerforationGravityDifference]) ≈ expected_gravity
    @test pu[:Well][:WellIndices] ≈ vcat(pc[:Well][:WellIndices], pc[:Well][:FractureWellIndices])
    @test !haskey(unified.models, :Fractures)

    # A hit at an existing well node is shared by its two adjacent segments.
    endpoint = setup_well(matrix, [1, 1, 2]; simple_well = false,
        well_cell_centers = [0.5 1 1.5; 0.5 0.5 0.5; 0.5 0.5 0.5])
    model = setup_embedded_fracture_model(matrix, domain, system;
        wells = [endpoint], extend_wells = true)
    w = physical_representation(model[:Well].data_domain)
    @test w.num_nodes == 3
    @test length(w.perforations.fracture) == 1
    @test only(w.perforations.self_fracture) == 2

    # A segment in a fracture plane crosses two clipped fracture cells.
    coplanar = setup_well(matrix, [1, 3]; simple_well = false, dir = :y,
        well_cell_centers = [1.0 1; 0.2 1.8; 0.5 0.5])
    model = setup_embedded_fracture_model(matrix, domain, system;
        wells = [coplanar], extend_wells = true)
    w = physical_representation(model[:Well].data_domain)
    @test length(w.perforations.fracture) == 2
    @test w.num_nodes == 4

    # Simultaneous hits on different fractures share a single well node.
    one = reservoir_domain(reservoir_mesh((1, 1, 1), (2.0, 2.0, 2.0)))
    objects = [fracture, RectangleFracture([1, 1, 1], [0, 1, 0], (2, 2))]
    crossing = setup_well(one, [1, 1]; simple_well = false,
        well_cell_centers = [0.2 1.8; 0.2 1.8; 0.5 0.5])
    model = setup_embedded_fracture_model(one, objects, system;
        wells = [crossing], extend_wells = true)
    w = physical_representation(model[:Well].data_domain)
    @test length(w.perforations.fracture) == 2
    @test w.num_nodes == 3
    @test w.perforations.self_fracture == [3, 3]

    miss = RectangleFracture([1, 1.8, 1.8], [1, 0, 0], (0.2, 0.2))
    model = setup_embedded_fracture_model(matrix, miss, system;
        wells = [well], extend_wells = true)
    @test !haskey(physical_representation(model[:Well].data_domain).perforations, :fracture)
end

@testset "EFM extended wells: coupled versus unified flow" begin
    for system in (ImmiscibleSystem((AqueousPhase(), LiquidPhase())), :geothermal)
        matrix = reservoir_domain(reservoir_mesh((3, 3, 3), (3.0, 3.0, 3.0));
            permeability = 1e-2si_unit(:darcy))
        fractures = [RectangleFracture([1, 1.5, 1.5], [1, 0, 0], (2.5, 2.5)),
            EllipticalFracture([2, 1.5, 1.5], [1, 0, 0], (1.2, 1.2))]
        inj = setup_well(matrix, [13, 15]; name = :Injector, simple_well = false, dir = :x, radius = 1e-2)
        prod = setup_well(matrix, [number_of_cells(matrix)]; name = :Producer, simple_well = false)
        wells = [inj, prod]
        models = [setup_fractured_reservoir_model(matrix, fractures, system;
            method = :efm, unified = unified, extend_wells = true, wells = wells)
            for unified in (false, true)]
        pc, pu = setup_parameters.(models)
        @test pu[:Injector][:WellIndices] ≈ vcat(pc[:Injector][:WellIndices], pc[:Injector][:FractureWellIndices])
        @test pu[:Injector][:PerforationGravityDifference] ≈ vcat(pc[:Injector][:PerforationGravityDifference], pc[:Injector][:FracturePerforationGravityDifference])
        @test length(physical_representation(models[1][:Injector].data_domain).perforations.fracture) == 2
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
            @test all(report -> length(report[:ministeps]) == 1 && only(report[:ministeps])[:success], result.result.reports)
            push!(results, result)
        end
        @test length(results[1].result.states) == length(results[2].result.states) == 4
        for (coupled, unified) in zip(results[1].result.states, results[2].result.states)
            @test unified[:Reservoir][:Pressure] ≈ vcat(coupled[:Reservoir][:Pressure], coupled[:Fractures][:Pressure]) rtol = 1e-7
            if key == :Temperature
                expected = vcat(coupled[:Reservoir][key], coupled[:Fractures][key])
            else
                expected = hcat(coupled[:Reservoir][key], coupled[:Fractures][key])
            end
            @test unified[:Reservoir][key] ≈ expected rtol = 1e-7 atol = 1e-7
        end
    end
end
