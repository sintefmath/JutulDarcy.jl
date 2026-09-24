using Jutul, JutulDarcy, Test

@testset "Merged wells" begin
    bar = si_unit(:bar)
    g = CartesianMesh((3, 3, 3), (300.0, 300.0, 60.0))
    domain = reservoir_domain(g, permeability = 0.1*si_unit(:darcy), porosity = 0.2)
    wells = [
        setup_vertical_well(domain, 1, 1, name = :SimpleProducer, simple_well = true),
        setup_vertical_well(domain, 3, 3, name = :SimpleInjector, simple_well = true),
        setup_vertical_well(domain, 1, 3, name = :SegmentProducer, simple_well = false),
        setup_vertical_well(domain, 3, 1, name = :SegmentInjector, simple_well = false),
    ]
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, sys, wells = wells, block_backend = false)
    state0 = setup_reservoir_state(model, Pressure = 200bar, Saturations = [0.8, 0.2])
    controls = Dict(
        :SimpleProducer => ProducerControl(BottomHolePressureTarget(180bar)),
        :SimpleInjector => InjectorControl(BottomHolePressureTarget(220bar), [1.0, 0.0], density = 1000.0),
        :SegmentProducer => ProducerControl(BottomHolePressureTarget(180bar)),
        :SegmentInjector => InjectorControl(BottomHolePressureTarget(220bar), [1.0, 0.0], density = 1000.0),
    )
    forces = setup_reservoir_forces(model, control = controls)
    forces[:SimpleProducer] = setup_forces(model[:SimpleProducer],
        mask = PerforationMask([1.0, 0.5, 1.0]))
    forces[:SegmentProducer] = setup_forces(model[:SegmentProducer],
        mask = PerforationMask([1.0, 0.5, 1.0]))
    case = JutulCase(model, [si_unit(:day)], forces, state0 = state0)
    merged = JutulDarcy.merge_similar_wells(case)
    @test Set(well_symbols(merged.model)) == Set(keys(controls))
    @test length(merged.model.models) == length(model.models) - 2
    @test physical_representation(merged.model.models[:SimpleWells]).multiwell.top_nodes == [1, 2]
    @test physical_representation(merged.model.models[:MultiSegmentWells]).multiwell.top_nodes == [1, 4]
    for source in (:SimpleWells, :MultiSegmentWells)
        multiwell = physical_representation(merged.model.models[source]).multiwell
        facility_names = merged.model.models[:Facility].domain.well_symbols
        expected_cells = [findfirst(isequal(name), facility_names) for name in multiwell.names]
        term_locations = (
            (JutulDarcy.WellFromFacilityFlowCT, source, :Facility),
            (JutulDarcy.FacilityFromWellBottomHolePressureCT, :Facility, source),
            (JutulDarcy.FacilityFromWellSurfaceComponentRatesCT, :Facility, source),
        )
        for (ct_type, target, source_model) in term_locations
            terms = [pair.cross_term for pair in merged.model.cross_terms
                if pair.target == target && pair.source == source_model &&
                   pair.cross_term isa ct_type]
            @test length(terms) == 1
            @test only(terms).wells == multiwell.names
            @test only(terms).facility_cells == expected_cells
            @test only(terms).well_cells == multiwell.top_nodes
        end
    end
    @test merged.state0[:SimpleWells][:Pressure] ==
        [state0[:SimpleProducer][:Pressure]; state0[:SimpleInjector][:Pressure]]
    well_domain = merged.model.models[:SimpleWells].data_domain
    well_model = SimulationModel(well_domain, JutulDarcy.SimpleWellSystem(sys))
    mass_equation = well_model.equations[:mass_conservation]
    convergence_storage = (
        state = (FluidVolume = [2.0, 4.0],),
        well_convergence_names = ["M1", "M2", "M3", "M4"],
    )
    residual = [1.0 2.0; -3.0 4.0]
    criterion = Jutul.convergence_criterion(well_model, convergence_storage,
        mass_equation, nothing, residual; dt = 2.0)
    @test criterion.CNV.errors ≈ [0.1, 0.3, 0.1, 0.2]
    @test residual == [1.0 2.0; -3.0 4.0]
    residual .= 0.0
    @test criterion.CNV.errors ≈ [0.1, 0.3, 0.1, 0.2]
    vector_storage = (
        state = convergence_storage.state,
        well_convergence_names = ["M1", "M2"],
    )
    vector_criterion = Jutul.convergence_criterion(well_model,
        vector_storage, mass_equation, nothing, [1.0, -2.0]; dt = 2.0)
    @test vector_criterion.CNV.errors ≈ [0.1, 0.1]
    fresh = setup_reservoir_state(merged.model, Pressure = 200bar,
        Saturations = [0.8, 0.2])
    @test fresh[:Facility][:BottomHolePressure] == state0[:Facility][:BottomHolePressure]

    ref = simulate_reservoir(case, info_level = -1)
    got = simulate_reservoir(merged, info_level = -1)
    got_cprw = simulate_reservoir(merged, precond = :cprw, info_level = -1)
    got_cprw_impes = simulate_reservoir(merged, precond = :cprw,
        linear_solver_arg = (cpr_type = :true_impes,), info_level = -1)
    for name in keys(controls)
        for (quantity, values) in ref.wells[name]
            other = got.wells[name][quantity]
            other_cprw = got_cprw.wells[name][quantity]
            other_cprw_impes = got_cprw_impes.wells[name][quantity]
            if eltype(values) <: Number
                @test isapprox(values, other; rtol = 1e-5, nans = true)
                @test isapprox(values, other_cprw; rtol = 1e-5, nans = true)
                @test isapprox(values, other_cprw_impes; rtol = 1e-5, nans = true)
            else
                @test values == other
                @test values == other_cprw
                @test values == other_cprw_impes
            end
        end
    end
    @test Set(keys(ref.summary["VALUES"]["WELLS"])) ==
        Set(keys(got.summary["VALUES"]["WELLS"]))
    for (name, outputs) in ref.summary["VALUES"]["WELLS"]
        for (quantity, values) in outputs
            @test isapprox(values, got.summary["VALUES"]["WELLS"][name][quantity];
                rtol = 1e-5, nans = true)
            @test isapprox(values, got_cprw.summary["VALUES"]["WELLS"][name][quantity];
                rtol = 1e-5, nans = true)
        end
    end
end

@testset "Merged CPRW with KernelAbstractions" begin
    bar = si_unit(:bar)
    g = CartesianMesh((2, 1, 2), (200.0, 100.0, 40.0))
    domain = reservoir_domain(g, permeability = 0.1*si_unit(:darcy), porosity = 0.2)
    wells = [
        setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = true),
        setup_vertical_well(domain, 2, 1, name = :Injector, simple_well = true),
    ]
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, sys, wells = wells, block_backend = true)
    state0 = setup_reservoir_state(model, Pressure = 200bar, Saturations = [0.8, 0.2])
    controls = Dict(
        :Producer => ProducerControl(BottomHolePressureTarget(190bar)),
        :Injector => InjectorControl(BottomHolePressureTarget(210bar),
            [1.0, 0.0], density = 1000.0),
    )
    forces = setup_reservoir_forces(model, control = controls)
    case = JutulDarcy.merge_similar_wells(JutulCase(model, [si_unit(:day)], forces, state0 = state0))
    cpu_model = setup_reservoir_model(domain, sys, wells = wells, block_backend = false)
    cpu_state0 = setup_reservoir_state(cpu_model, Pressure = 200bar,
        Saturations = [0.8, 0.2])
    cpu_forces = setup_reservoir_forces(cpu_model, control = controls)
    cpu_case = JutulDarcy.merge_similar_wells(JutulCase(cpu_model, [si_unit(:day)],
        cpu_forces, state0 = cpu_state0))
    cpu = simulate_reservoir(cpu_case, precond = :cprw, info_level = -1)
    ka = simulate_reservoir(case, mode = :ka, precond = :cprw, info_level = -1)
    for name in keys(controls)
        for quantity in (:bhp, :mass_rate)
            @test isapprox(cpu.wells[name][quantity], ka.wells[name][quantity]; rtol = 1e-5)
        end
    end
    add_tracers_to_model!(case.model, SinglePhaseTracer(1))
    tracer_terms = [pair.cross_term for pair in case.model.cross_terms
        if pair.target == :SimpleWells && pair.source == :Facility &&
           pair.cross_term isa JutulDarcy.Tracers.WellFromFacilityTracerCT]
    @test length(tracer_terms) == 1
    @test only(tracer_terms).wells == [:Producer, :Injector]
    @test only(tracer_terms).well_cells ==
        physical_representation(case.model.models[:SimpleWells]).multiwell.top_nodes
end

@testset "Merged thermal wells" begin
    bar = si_unit(:bar)
    g = CartesianMesh((2, 2, 2), (200.0, 200.0, 40.0))
    domain = reservoir_domain(g, permeability = 0.1*si_unit(:darcy), porosity = 0.2)
    wells = [
        setup_vertical_well(domain, 1, 1, name = :ThermalProducer, simple_well = true),
        setup_vertical_well(domain, 2, 2, name = :ThermalInjector, simple_well = true),
    ]
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, sys, wells = wells, thermal = true,
        block_backend = false)
    state0 = setup_reservoir_state(model, Pressure = 200bar,
        Saturations = [0.8, 0.2], Temperature = 300.0)
    controls = Dict(
        :ThermalProducer => ProducerControl(BottomHolePressureTarget(190bar)),
        :ThermalInjector => InjectorControl(BottomHolePressureTarget(210bar),
            [1.0, 0.0], density = 1000.0, temperature = 320.0),
    )
    forces = setup_reservoir_forces(model, control = controls)
    case = JutulCase(model, [si_unit(:day)], forces, state0 = state0)
    merged = JutulDarcy.merge_similar_wells(case)
    multiwell = physical_representation(merged.model.models[:SimpleWells]).multiwell
    term_locations = (
        (JutulDarcy.WellFromFacilityThermalCT, :SimpleWells, :Facility),
        (JutulDarcy.FacilityFromWellTemperatureCT, :Facility, :SimpleWells),
        (JutulDarcy.FacilityFromWellEnthalpyCT, :Facility, :SimpleWells),
    )
    for (ct_type, target, source_model) in term_locations
        terms = [pair.cross_term for pair in merged.model.cross_terms
            if pair.target == target && pair.source == source_model &&
               pair.cross_term isa ct_type]
        @test length(terms) == 1
        @test only(terms).wells == multiwell.names
        @test only(terms).well_cells == multiwell.top_nodes
    end
    ref = simulate_reservoir(case, info_level = -1)
    got = simulate_reservoir(merged, info_level = -1)
    got_cprw = simulate_reservoir(merged, precond = :cprw, info_level = -1)
    for name in keys(controls)
        for quantity in (:bhp, :mass_rate, :temperature)
            @test isapprox(ref.wells[name][quantity], got.wells[name][quantity]; rtol = 1e-5)
            @test isapprox(ref.wells[name][quantity], got_cprw.wells[name][quantity]; rtol = 1e-5)
        end
    end
end
