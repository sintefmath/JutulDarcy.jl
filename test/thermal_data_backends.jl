using Test, Jutul, JutulDarcy, Adapt, CUDA
using JutulDarcy.StaticArrays: SVector
using Jutul.KernelExecution: KernelAbstractions
using .KernelAbstractions: @kernel, @index

@kernel function injection_enthalpy_test_kernel!(out, evaluator, model, state, mixture, temperatures)
    cell = @index(Global)
    ctrl = (injection_mixture = mixture,)
    out[cell] = JutulDarcy.injected_fluid_enthalpy(ctrl, evaluator, model,
        state, state.Pressure[cell], temperatures[cell], cell)
end

@testset "Injection enthalpy evaluates without model definitions" begin
    if CUDA.functional()
        CUDA.allowscalar(false)
    end
    system = ImmiscibleSystem((AqueousPhase(), VaporPhase()))
    model = (system = system,)
    temperatures = [290.0, 320.0]
    state = (Pressure = [1e5, 2e5], PhaseMassDensities = [1000.0 990.0; 1.0 2.0],
        ComponentHeatCapacity = [1500.0 1500.0; 500.0 500.0])
    mixture = SVector(0.8, 0.2)
    capacities = [SVector(1000.0, 2000.0), SVector(2000.0, 3000.0)]
    energy = JutulDarcy.IntegratedHeatCapacity([280.0, 340.0], capacities)
    variable = JutulDarcy.TemperatureDependentVariable((energy, energy); regions = [1, 2])
    for pressure_work in (Val(false), Val(true))
        evaluator = JutulDarcy.InjectionEnthalpyEvaluator(variable, pressure_work)
        # The mass-weighted heat capacity is 1200 at 280 K and increases
        # linearly with slope 1000/60. Integrate it analytically.
        reference = [336000 + 1200*dT + (1000/120)*dT^2 for dT in (10.0, 40.0)]
        if pressure_work === Val(true)
            reference .+= [1e5*(0.8/1000 + 0.2), 2e5*(0.8/990 + 0.2/2)]
        end
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
end

@testset "Thermal DATA cases on ka_cuda" begin
    flush(stdout)
    if CUDA.functional()
        CUDA.allowscalar(false)
        fixtures = joinpath(@__DIR__, "data", "thermal")
        for name in ("HEAT_SMALL.DATA", "THERMAL_TRACER_SMALL.DATA", "TRACER_SMALL.DATA")
            @testset "$name" begin
                case = setup_case_from_data_file(joinpath(fixtures, name);
                    dp_max_abs_well = Inf, parse_arg = (silent = true,))
                options = (info_level = -1, error_on_incomplete = true,
                    nonlinear_tolerance = 1e-7, tol_cnv = 1e-6, tol_mb = 1e-9,
                    tol_cnv_well = 1e-7)
                reference = simulate_reservoir(case; options...)
                # The default cuSPARSE block ILU does not pivot. This dry-water
                # startup needs a cell-block method that handles zero scalar
                # pivots; use the existing portable block DILU smoother.
                result = simulate_reservoir(case; mode = :ka_cuda,
                    linear_solver_arg = (smoother_type = :ka_dilu,), options...)
                @test length(result.states) == length(case.dt)
                for (actual, expected) in zip(result.states, reference.states)
                    @test actual[:Pressure] ≈ expected[:Pressure] rtol = 1e-5
                    for property in (:Temperature, :TracerConcentrations)
                        if haskey(expected, property)
                            @test actual[property] ≈ expected[property] rtol = 1e-4 atol = 1e-7
                        end
                    end
                end
            end
        end
    else
        @test_skip CUDA.functional()
    end
end
