using Jutul, JutulDarcy, GeoEnergyIO, Test, CUDA

@testset "Tracer injection follows the prescribed carrier stream" begin
    ctrl = InjectorControl(TotalRateTarget(1.0), [1.0, 0.0, 0.0]; tracers = [0.8])
    cfg = JutulDarcy.WellGroupConfiguration([:I], Dict(:I => ctrl))
    facility_state = (WellGroupConfiguration = cfg, TotalSurfaceMassRate = [2.0])
    # The well initially contains gas, but the injected stream is pure water.
    well_state = (Saturations = reshape([0.0, 0.0, 1.0], 3, 1),
        PhaseMassDensities = reshape([1000.0, 800.0, 1.0], 3, 1))
    out = zeros(1)
    for system in (ImmiscibleSystem(:wog), StandardBlackOilSystem())
        flux = JutulDarcy.Tracers.TracerFluxType(SinglePhaseTracer(system, AqueousPhase()))
        ct = JutulDarcy.Tracers.WellFromFacilityTracerCT([:I], [1], [1], flux.tracers)
        Jutul.update_cross_term_in_entity!(out, 1, well_state, nothing,
            facility_state, nothing, (system = system,), nothing, ct, nothing, 1.0, nothing)
        @test only(out) ≈ -1.6
    end
end

@testset "Thermal and passive tracer DATA regressions" begin
    fixtures = joinpath(@__DIR__, "data", "thermal")
    # Use the same timestep sequence for CPU and backend comparisons.
    options = (info_level = -1, error_on_incomplete = true, timesteps = :none,
        nonlinear_tolerance = 1e-7, tol_cnv = 1e-6, tol_mb = 1e-9,
        tol_cnv_well = 1e-7)
    if CUDA.functional()
        CUDA.allowscalar(false)
    end
    for (name, nc, thermal, tracer) in [
            ("HEAT_SMALL.DATA", 50, true, false),
            ("THERMAL_TRACER_SMALL.DATA", 105, true, true),
            ("TRACER_SMALL.DATA", 25, false, true)
        ]
        @testset "$name" begin
            case = setup_case_from_data_file(joinpath(fixtures, name);
                parse_arg = (silent = true,), dp_max_abs_well = Inf)
            rmodel = reservoir_model(case.model)
            @test number_of_cells(rmodel.data_domain) == nc
            @test number_of_tracers(case.model) == Int(tracer)
            if thermal
                @test Jutul.minimum_value(rmodel[:Temperature]) == 0.0
            end
            if thermal && !tracer
                units = GeoEnergyIO.InputParser.DeckUnitSystem(:field)
                properties = Jutul.evaluate_all_secondary_variables(rmodel,
                    case.state0[:Reservoir], case.parameters[:Reservoir])
                @test properties[:FluidEnthalpy] ≈ properties[:FluidInternalEnergy] .+
                    properties[:Pressure]' ./ properties[:PhaseMassDensities]
                @test vec(properties[:RockHeatCapacity]) .* case.parameters[:Reservoir][:RockDensity] ≈ fill(35*units.volume_heat_capacity, nc)
                @test rmodel.data_domain[:bulk_thermal_conductivity] ≈ fill(24*units.rock_conductivity, nc)
                @test case.forces[1][:Facility].control[:I1].temperature ≈ convert_to_si(20.0, :Fahrenheit)
            end
            if tracer
                @test all(iszero, case.state0[:Reservoir][:TracerConcentrations])
                if thermal
                    injector = :I1
                else
                    injector = :I
                end
                @test case.forces[1][:Facility].control[injector].tracers == [1.0]
                if thermal
                    ctrl = case.forces[end][:Facility].control[injector]
                    @test ctrl.temperature ≈ 363.15
                    @test ctrl.tracers == [0.5]
                end
            end
            result = simulate_reservoir(case; options...)
            @test length(result.states) == length(case.dt)
            for state in result.states
                @test all(p -> isfinite(p) && p > 0, state[:Pressure])
                if thermal
                    @test all(T -> isfinite(T) && T > 0, state[:Temperature])
                end
                if tracer
                    concentrations = state[:TracerConcentrations]
                    @test all(C -> isfinite(C) && C >= -1e-8, concentrations)
                    # Concentration is ill-conditioned where the carrier phase
                    # has vanished. Check the phase-weighted tracer amount,
                    # using the mass-balance tolerance rather than dry-cell C.
                    if haskey(state, :Saturations)
                        carrier = state[:Saturations][1, :]
                        @test maximum((vec(concentrations) .- 1).*carrier) <= 1e-7
                    else
                        @test maximum(concentrations) <= 1 + 1e-7
                    end
                end
            end
            if tracer
                @test maximum(result.states[end][:TracerConcentrations]) > 0
                if thermal
                    @test case.forces[1][:Facility].control[:P1].target isa SurfaceWaterRateTarget
                    @test first(result.wells[:P1, :bhp]) ≈ si_unit(:bar) rtol = 1e-6
                end
            end
            if thermal
                @test minimum(result.states[end][:Temperature]) < minimum(case.state0[:Reservoir][:Temperature])
            end
            @testset "ka_cuda" begin
                if CUDA.functional()
                    policies = NamedTuple[NamedTuple()]
                    if thermal && tracer
                        push!(policies, (wells_on_device = true, mixed_cross_terms_on_host = true))
                    end
                    for policy in policies
                        # Block DILU handles the zero scalar pivots at dry-water startup.
                        gpu = simulate_reservoir(case; mode = :ka_cuda, policy...,
                            linear_solver_arg = (smoother_type = :ka_dilu,), options...)
                        @test length(gpu.states) == length(case.dt)
                        for (actual, expected) in zip(gpu.states, result.states)
                            @test actual[:Pressure] ≈ expected[:Pressure] rtol = 1e-5
                            for (property, tolerance) in ((:Temperature, 1e-5), (:TracerMasses, 1e-4))
                                if haskey(expected, property)
                                    @test actual[property] ≈ expected[property] rtol = tolerance atol = 1e-7
                                end
                            end
                        end
                    end
                else
                    @test_skip CUDA.functional()
                end
            end
        end
    end
end
