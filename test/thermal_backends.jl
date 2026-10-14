using Test, Jutul, JutulDarcy, Adapt, JLArrays

@testset "Tabulated thermal properties on JLArrays" begin
    for physics in (:geothermal, :co2brine_immiscible, :co2brine_compositional)
        @testset "$physics" begin
            domain = get_1d_reservoir(3)
            is_geothermal = physics == :geothermal
            is_compositional = physics == :co2brine_compositional
            options = is_geothermal ? NamedTuple() :
                (co2_physics = is_compositional ? :kvalue : :immiscible,
                    co2_source = :csp11)
            model, parameters = setup_reservoir_model(domain,
                is_geothermal ? :geothermal : :co2brine;
                thermal = true, extra_out = true, block_backend = false, options...)
            composition = is_geothermal ? NamedTuple() : is_compositional ?
                (OverallMoleFractions = [0.8, 0.2],) : (Saturations = [0.8, 0.2],)
            state0 = setup_reservoir_state(model;
                Pressure = [2.1e7, 2.0e7, 1.9e7],
                Temperature = [330.0, 325.0, 320.0], composition...)
            model = reservoir_model(model)
            parameters = parameters[:Reservoir]
            state0 = state0[:Reservoir]
            reference_simulator = Simulator(model; state0, parameters)
            reference, = simulate!(reference_simulator, [1.0]; info_level = -1)
            simulator = transfer_to_backend(Simulator(model; state0, parameters), JLBackend())
            enthalpy = simulator.model.secondary_variables.FluidEnthalpy
            @test enthalpy.tab[1].F isa JLArray
            @test simulator.model.secondary_variables.FluidInternalEnergy isa
                JutulDarcy.FluidInternalEnergyFromFluidEnthalpy
            states, = simulate!(simulator, [1.0]; info_level = -1)
            @test length(states) == length(reference) == 1
            for name in (:Pressure, :Temperature, :FluidEnthalpy, :FluidInternalEnergy)
                @test Array(states[end][name]) ≈ reference[end][name] rtol = 1e-8
            end
        end
    end

    # Backend transfer must retain both regional lookup and the table length.
    tables = JutulDarcy.CO2Properties.co2_brine_property_tables(co2_source = :csp11)
    enthalpy = JutulDarcy.PressureTemperatureDependentEnthalpy(
        (tables[:enthalpy], tables[:enthalpy]); regions = [1, 2, 1])
    context = KernelAbstractionsContext(JLBackend(); float_type = Float32, index_type = Int32)
    adapted = Adapt.adapt(context, enthalpy)
    @test adapted.regions isa JLArray{Int32}
    @test Array(adapted.regions) == [1, 2, 1]
    @test eltype(eltype(adapted.tab[1].F)) === Float32
end

@testset "Tabulated geothermal wells on a KA backend" begin
    case = JutulDarcy.setup_mini_wellcase(Val(:geothermal);
        simple_well = false, nstep = 1, total_time = 0.01*si_unit(:day))
    for (backend, float_type) in (
            (JLBackend(), Float64),
            (Jutul.KernelExecution.KernelAbstractions.CPU(), Float32))
        simulator, = setup_reservoir_simulator(case;
            mode = :ka, ka_backend = backend, float_type,
            info_level = -1, linear_solver = nothing, timesteps = :none)
        forces = Jutul.preprocess_forces(simulator, case.forces).forces
        dt = only(case.dt)
        Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
        Jutul.update_state_dependents!(
            simulator.storage, simulator.model, dt, forces; time = dt)
        Jutul.update_linearized_system!(simulator.storage, simulator.model)
        for name in (:Reservoir, :Injector, :Producer)
            state = simulator.storage[name].state
            @test all(isfinite, Jutul.value.(Array(state.FluidEnthalpy)))
            @test all(isfinite, Jutul.value.(Array(state.FluidInternalEnergy)))
        end
        @test all(isfinite, Array(simulator.storage.LinearizedSystem.r_buffer))
    end
end
