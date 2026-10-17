using Jutul, JutulDarcy, GeoEnergyIO, Test

@testset "Thermal well enthalpy uses mass fractions" begin
    state = (Saturations = reshape([0.5, 0.5], 2, 1),
        PhaseMassDensities = reshape([1000.0, 1.0], 2, 1),
        FluidEnthalpy = reshape([1e5, 1e6], 2, 1),
        ComponentHeatCapacity = reshape([1500.0, 500.0], 2, 1))
    @test JutulDarcy.well_top_node_enthalpy(nothing, state, 1) ≈ (1000e5 + 1e6)/1001
    system = ImmiscibleSystem((AqueousPhase(), VaporPhase()))
    ctrl = InjectorControl(TotalRateTarget(1.0), [1.0, 0.0])
    evaluator = JutulDarcy.InjectionEnthalpyEvaluator(nothing, Val(true))
    H = JutulDarcy.injected_fluid_enthalpy(ctrl, evaluator,
        (system = system,), state, 1e5, 300.0, 1)
    @test H ≈ 1500*300 + 1e5/1000
end

@testset "Integrate temperature-dependent heat capacities" begin
    table = JutulDarcy.IntegratedHeatCapacity([280.0, 300.0, 340.0], [1000.0, 2000.0, 3000.0])
    @test table(300.0) ≈ 310000.0
    @test table(320.0) ≈ 355000.0
    for (T, capacity) in [(270.0, 1000.0), (290.0, 1500.0), (320.0, 2500.0), (350.0, 3000.0)]
        @test Jutul.ForwardDiff.derivative(table, T) ≈ capacity
    end
end

@testset "Tracer injection follows the prescribed carrier stream" begin
    system = ImmiscibleSystem((AqueousPhase(), VaporPhase()))
    flux = JutulDarcy.Tracers.TracerFluxType(SinglePhaseTracer(system, AqueousPhase()))
    well = (system = system,)
    ctrl = InjectorControl(TotalRateTarget(1.0), [1.0, 0.0]; tracers = [0.8])
    cfg = JutulDarcy.WellGroupConfiguration([:I], Dict(:I => ctrl))
    facility_state = (WellGroupConfiguration = cfg, TotalSurfaceMassRate = [2.0])
    # The well initially contains gas, but the injected stream is pure water.
    well_state = (Saturations = reshape([0.0, 1.0], 2, 1),
        PhaseMassDensities = reshape([1000.0, 1.0], 2, 1))
    ct = JutulDarcy.Tracers.WellFromFacilityTracerCT([:I], [1], [1], flux.tracers)
    out = zeros(1)
    Jutul.update_cross_term_in_entity!(out, 1, well_state, nothing,
        facility_state, nothing, well, nothing, ct, nothing, 1.0, nothing)
    @test only(out) ≈ -1.6
end

@testset "Thermal and passive tracer DATA regressions" begin
    fixtures = joinpath(@__DIR__, "data", "thermal")
    for (name, nc, thermal, tracer) in [
            ("HEAT_SMALL.DATA", 50, true, false),
            ("THERMAL_TRACER_SMALL.DATA", 105, true, true),
            ("TRACER_SMALL.DATA", 25, false, true)
        ]
        @testset "$name" begin
            case, data = setup_case_from_data_file(joinpath(fixtures, name);
                include_data = true, parse_arg = (silent = true,), dp_max_abs_well = Inf)
            rmodel = reservoir_model(case.model)
            @test number_of_cells(rmodel.data_domain) == nc
            @test number_of_tracers(case.model) == Int(tracer)
            if thermal
                @test all(isfinite, case.state0[:Reservoir][:Temperature])
                @test all(isfinite, case.state0[:Facility][:SurfaceEnthalpy])
                properties = Jutul.evaluate_all_secondary_variables(rmodel,
                    case.state0[:Reservoir], case.parameters[:Reservoir])
                @test all(properties[:ComponentHeatCapacity] .> 0)
                @test all(properties[:RockHeatCapacity] .> 0)
                @test all(properties[:PhaseViscosities] .> 0)
                @test properties[:FluidEnthalpy] == properties[:FluidInternalEnergy]
                if !tracer
                    @test :Temperature in Jutul.get_dependencies(rmodel.secondary_variables[:PhaseViscosities], rmodel)
                    @test vec(properties[:RockHeatCapacity]) .* case.parameters[:Reservoir][:RockDensity] ≈ fill(35000.0, nc)
                    @test rmodel.data_domain[:rock_thermal_conductivity] ≈ fill(24/86.4, nc)
                    @test case.parameters[:Reservoir][:RockThermalConductivities] ≈ JutulDarcy.compute_face_trans(rmodel.data_domain, fill(24/86.4, nc))
                end
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
            result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true,
                nonlinear_tolerance = 1e-7, tol_cnv = 1e-6, tol_mb = 1e-9,
                tol_cnv_well = 1e-7)
            @test length(result.states) == length(case.dt)
            for state in result.states
                @test all(isfinite, state[:Pressure])
                @test all(state[:Pressure] .> 0)
                if thermal
                    @test all(isfinite, state[:Temperature])
                    @test minimum(state[:Temperature]) > 273.15
                end
                if tracer
                    concentrations = state[:TracerConcentrations]
                    @test all(isfinite, concentrations)
                    @test minimum(concentrations) >= -1e-8
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
            end
            if thermal
                @test minimum(result.states[end][:Temperature]) < minimum(case.state0[:Reservoir][:Temperature])
            end
        end
    end
end
