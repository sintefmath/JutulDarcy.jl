import JutulDarcy.EFM

@testset "Unified DFM domains and coefficients" begin
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
    for strategy in (:star_delta, :keep, :remove)
        d = setup_cross_dfm(intersection_strategy = strategy)
        before_matrix = deepcopy(d.matrix.data)
        before_fractures = deepcopy(d.fractures.data)
        coupled = setup_fractured_reservoir_model(EFM.copy_fracture_data(d.matrix),
            EFM.copy_fracture_data(d.fractures), system; thermal = true)
        pc = setup_parameters(coupled)
        repeated = setup_parameters(coupled)
        @test repeated[:Fractures][:FractureMatrixThermalConductivity] == pc[:Fractures][:FractureMatrixThermalConductivity]
        @test repeated[:Fractures][:RockThermalConductivities] == pc[:Fractures][:RockThermalConductivities]
        @test repeated[:Fractures][:FluidThermalConductivities] == pc[:Fractures][:FluidThermalConductivities]
        unified, pu = setup_fractured_reservoir_model(d.matrix, d.fractures, system;
            unified = true, thermal = true, extra_out = true)
        @test !haskey(unified.models, :Fractures)
        @test !haskey(pu, :Fractures)
        domain = reservoir_domain(unified)
        nm = number_of_cells(d.matrix)
        nf = number_of_cells(d.fractures)
        @test number_of_cells(domain) == nm+nf
        @test domain[:is_fracture] isa Vector{Int}
        @test domain[:is_fracture] == vcat(zeros(Int, nm), ones(Int, nf))
        @test domain[:volumes] ≈ vcat(reservoir_domain(coupled)[:volumes], coupled[:Fractures].data_domain[:volumes])
        @test domain[:porosity] == vcat(d.matrix[:porosity], d.fractures[:porosity])
        @test domain[:cell_centroids] == hcat(d.matrix[:cell_centroids], d.fractures[:cell_centroids])
        @test pu[:Reservoir][:Transmissibilities] ≈ vcat(pc[:Reservoir][:Transmissibilities],
            pc[:Fractures][:Transmissibilities], pc[:Fractures][:FractureMatrixTransmissibility])
        @test pu[:Reservoir][:RockThermalConductivities] ≈ vcat(pc[:Reservoir][:RockThermalConductivities],
            pc[:Fractures][:RockThermalConductivities], pc[:Fractures][:FractureMatrixThermalConductivity])
        @test pu[:Reservoir][:FluidThermalConductivities] ≈ hcat(pc[:Reservoir][:FluidThermalConductivities],
            pc[:Fractures][:FluidThermalConductivities], zeros(2, length(d.fractures[:matrix_cells])))
        # Setup must preserve input geometry and properties, including heat data.
        @test isequal(d.matrix.data, before_matrix)
        @test isequal(d.fractures.data, before_fractures)
        unadjusted = setup_fractured_reservoir_model(d.matrix, d.fractures, system;
            unified = true, subtract_fracture_volume = false)
        @test reservoir_domain(unadjusted)[:volumes][1:nm] == d.matrix[:volumes]
        mesh = physical_representation(domain)
        @test mesh isa EmbeddedReservoirMesh
        @test occursin("DFM", sprint(show, MIME"text/plain"(), mesh))
        @test length(sprint(show, mesh)) < 1000
        converted = UnstructuredMesh(mesh)
        @test converted isa UnstructuredMesh{3}
        @test number_of_cells(converted) == nm+nf
        @test get_neighborship(converted) == domain[:neighbors]
        for cell in nm+1:nm+nf
            center, volume = Jutul.compute_centroid_and_measure(converted, Cells(), cell)
            @test volume ≈ domain[:volumes][cell] rtol = 1e-8
            @test center ≈ domain[:cell_centroids][:, cell] atol = 1e-8
            @test cell_ijk(mesh, cell) == cell_ijk(converted, cell)
        end
        for flatten in (true, false)
            plotted = Jutul.triangulate_mesh(mesh; flatten = flatten)
            @test maximum(plotted.mapper.indices.Cells) == nm+nf
        end
    end
end

@testset "Unified DFM flow and fracture well completions" begin
    for system in (SinglePhaseSystem(AqueousPhase()),
            ImmiscibleSystem((AqueousPhase(), LiquidPhase())), :geothermal)
        thermal = system == :geothermal
        for extend_wells in (false, true)
            d = setup_cross_dfm(intersection_strategy = :keep)
            # A conductive permeability with moderate contrast keeps the tiny
            # intersection cells above the double-precision residual floor.
            d.fractures[:permeability] .= 100si_unit(:darcy)
            # This well crosses a fracture between adjacent matrix cells.
            cells = d.matrix[:neighbors][:, first(d.fracture_mesh.parent_faces)]
            axis = argmax(abs.(d.matrix[:normals][:, first(d.fracture_mesh.parent_faces)]))
            direction = (:x, :y, :z)[axis]
            inj = setup_well(d.matrix, cells; name = :Injector, simple_well = false,
                dir = direction, radius = 1e-3)
            prod = setup_well(d.matrix, [number_of_cells(d.matrix)]; name = :Producer, simple_well = true, radius = 1e-3)
            wells = [inj, prod]
            before = deepcopy(inj.data)
            coupled = setup_fractured_reservoir_model(EFM.copy_fracture_data(d.matrix),
                EFM.copy_fracture_data(d.fractures), system;
                wells = deepcopy(wells), extend_wells = extend_wells)
            unified, pu = setup_fractured_reservoir_model(d.matrix, d.fractures, system;
                unified = true, wells = wells, extend_wells = extend_wells, extra_out = true)
            @test isequal(inj.data, before)
            pc = setup_parameters(coupled)
            cw = physical_representation(coupled[:Injector].data_domain)
            uw = physical_representation(unified[:Injector].data_domain)
            nm = number_of_cells(d.matrix)
            if extend_wells
                @test !isempty(cw.perforations.fracture)
                @test uw.perforations.reservoir == vcat(cw.perforations.reservoir, cw.perforations.fracture .+ nm)
                @test pu[:Injector][:WellIndices] ≈ vcat(pc[:Injector][:WellIndices], pc[:Injector][:FractureWellIndices])
                if thermal
                    @test pu[:Injector][:WellIndicesThermal] ≈ vcat(pc[:Injector][:WellIndicesThermal], pc[:Injector][:FractureWellIndicesThermal])
                end
            else
                @test uw.perforations.reservoir == cw.perforations.reservoir
            end
            initial = (Pressure = 50bar,)
            injection = (density = 1000.0,)
            mix = [1.0]
            keys_to_compare = [:Pressure]
            if system isa ImmiscibleSystem && number_of_phases(system) == 2
                initial = (initial..., Saturations = [0.0, 1.0])
                mix = [1.0, 0.0]
                push!(keys_to_compare, :Saturations)
            end
            if thermal
                initial = (initial..., Temperature = 363.15)
                injection = (injection..., temperature = 283.15)
                push!(keys_to_compare, :Temperature)
            end
            controls = Dict(:Injector => InjectorControl(TotalRateTarget(1e-5), mix; injection...),
                :Producer => ProducerControl(BottomHolePressureTarget(10bar)))
            function run(model)
                state0 = setup_reservoir_state(model; initial...)
                forces = setup_reservoir_forces(model; control = controls)
                case = JutulCase(model, fill(20.0, 4), forces; state0 = state0)
                solver = GenericKrylov(preconditioner = ILUZeroPreconditioner(),
                    relative_tolerance = 1e-10, absolute_tolerance = 1e-16)
                result = simulate_reservoir(case; info_level = -1, max_nonlinear_iterations = 40,
                    linear_solver = solver, tol_cnv = 1e-8, tol_mb = 1e-9, tol_cnve = 1e-8,
                    tol_cnv_well = 1e-6, tol_cnve_well = 1e-6,
                    timesteps = :none, max_timestep_cuts = 0)
                @test all(report -> length(report[:ministeps]) == 1 && only(report[:ministeps])[:success], result.result.reports)
                return result
            end
            rc, ru = run(coupled), run(unified)
            @test length(rc.result.states) == length(ru.result.states) == 4
            for (sc, su) in zip(rc.result.states, ru.result.states)
                for key in keys_to_compare
                    left, right = sc[:Reservoir][key], sc[:Fractures][key]
                    if left isa AbstractMatrix
                        expected = hcat(left, right)
                    else
                        expected = vcat(left, right)
                    end
                    @test su[:Reservoir][key] ≈ expected rtol = 1e-7 atol = 1e-7
                end
            end
            @test maximum(abs.(ru.result.states[end][:Reservoir][:Pressure] .- 50bar)) > 1.0
        end
    end
end
