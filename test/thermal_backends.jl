using Test, Jutul, JutulDarcy, Adapt, CUDA
using JutulDarcy.StaticArrays: SVector
using Jutul.KernelExecution: KernelAbstractions
using .KernelAbstractions: @kernel, @index

CUDA.allowscalar(false)

@testset "Tabulated thermal properties on POCL" begin
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
            simulator = transfer_to_backend(Simulator(model; state0, parameters), KernelAbstractions.CPU())
            enthalpy = simulator.model.secondary_variables.FluidEnthalpy
            @test enthalpy.tab[1].F isa Array
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
    context = KernelAbstractionsContext(KernelAbstractions.CPU(); float_type = Float32, index_type = Int32)
    adapted = Adapt.adapt(context, enthalpy)
    @test adapted.regions isa Array{Int32}
    @test Array(adapted.regions) == [1, 2, 1]
    @test eltype(eltype(adapted.tab[1].F)) === Float32
end

@testset "Tabulated geothermal wells on a KA backend" begin
    case = JutulDarcy.setup_mini_wellcase(Val(:geothermal);
        simple_well = false, nstep = 1, total_time = 0.01*si_unit(:day))
    case = JutulCase(convert_to_immutable_storage(case.model), case.dt, case.forces;
        state0 = case.state0, parameters = case.parameters)
    for (backend, float_type) in (
            (KernelAbstractions.CPU(), Float64),
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

@kernel function injection_enthalpy_test_kernel!(out, evaluator, model, state, mixture, temperatures)
    cell = @index(Global)
    ctrl = (injection_mixture = mixture,)
    out[cell] = JutulDarcy.injected_fluid_enthalpy(ctrl, evaluator, model,
        state, state.Pressure[cell], temperatures[cell], cell)
end

@testset "Injection enthalpy evaluates without model definitions" begin
    system = ImmiscibleSystem((AqueousPhase(), VaporPhase()))
    model = (system = system,)
    temperatures = [290.0, 320.0]
    state = (Pressure = [1e5, 2e5], PhaseMassDensities = [1000.0 990.0; 1.0 2.0],
        ComponentHeatCapacity = [1500.0 1500.0; 500.0 500.0])
    mixture = SVector(0.8, 0.2)
    capacities = [SVector(1000.0, 2000.0), SVector(2000.0, 3000.0)]
    energy = JutulDarcy.IntegratedHeatCapacity([280.0, 340.0], capacities)
    variable = JutulDarcy.TemperatureDependentVariable((energy, energy); regions = [1, 2])
    evaluator = JutulDarcy.InjectionEnthalpyEvaluator(variable)
    # The mass-weighted heat capacity is 1200 at 280 K and increases
    # linearly with slope 1000/60. Integrate it analytically and add p/ρ.
    reference = [336000 + 1200*dT + (1000/120)*dT^2 for dT in (10.0, 40.0)]
    reference .+= [1e5*(0.8/1000 + 0.2), 2e5*(0.8/990 + 0.2/2)]
    backends = Any[KernelAbstractions.CPU()]
    if CUDA.functional()
        push!(backends, CUDA.CUDABackend(always_inline = true))
    end
    for backend in backends, float_type in (Float32, Float64)
        context = KernelAbstractionsContext(backend; float_type)
        device_state = Adapt.adapt(context, state)
        device_evaluator = Adapt.adapt(context, evaluator)
        output = Adapt.adapt(context, zeros(2))
        launch! = injection_enthalpy_test_kernel!(backend, 32)
        launch!(output, device_evaluator, model, device_state, Adapt.adapt(context, mixture),
            Adapt.adapt(context, temperatures); ndrange = 2)
        KernelAbstractions.synchronize(backend)
        @test Array(output) ≈ reference rtol = 1e-6
    end
end
