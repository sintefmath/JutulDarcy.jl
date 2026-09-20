using Jutul
using JutulDarcy
using Test

const f32_backend = Jutul.KernelExecution.KernelAbstractions.CPU()
const f32_step = 0.01*si_unit(:day)

function short_reservoir_case(case)
    return JutulCase(case.model, [f32_step], case.forces[1:1];
        state0 = case.state0, parameters = case.parameters)
end

function test_float32_reservoir_case(case, linear_float, linear_index;
        precond = :cpr)
    simulator, config = setup_reservoir_simulator(case;
        mode = :ka,
        ka_backend = f32_backend,
        float_type = Float32,
        index_type = Int32,
        linear_float_type = linear_float,
        linear_index_type = linear_index,
        precond = precond,
        info_level = -1,
        failure_cuts_timestep = false,
        error_on_incomplete = true)

    context = simulator.model.context
    system = simulator.storage.LinearizedSystem
    pressure = simulator.storage.Reservoir.state.Pressure
    @test Jutul.float_type(context) === Float32
    @test Jutul.index_type(context) === Int32
    @test Jutul.linear_float_type(context) === linear_float
    @test Jutul.linear_index_type(context) === linear_index
    @test typeof(Jutul.value(first(pressure))) === Float32
    @test eltype(system.r_buffer) === linear_float
    @test eltype(system[1, 1].jac.rowptr) === linear_index
    @test config[:tolerances][:Facility][:control_equation].Abs == 0.01

    result = simulate_reservoir(case; simulator, config)
    solver = config[:linear_solver]
    @test eltype(solver.storage.x) === linear_float
    if solver.preconditioner isa CPRPreconditioner
        @test eltype(solver.preconditioner.storage.w_rhs) === linear_float
    end
    @test length(result.states) == 1
    if length(result.states) == 1
        @test all(isfinite, result.states[1][:Pressure])
    end
end

@testset "Float32 reservoir deck steps" begin
    for name in ("EGG", "SPE1", "SPE9")
        @testset "$name" begin
            path = JutulDarcy.GeoEnergyIO.test_input_file_path(
                name, "$name.DATA")
            case = setup_case_from_data_file(path;
                backend = :csr, block_backend = true)[1:1]
            short_case = short_reservoir_case(case)
            for (linear_float, linear_index) in
                    ((Float32, Int32), (Float64, Int64))
                @testset "linear $linear_float/$linear_index" begin
                    test_float32_reservoir_case(
                        short_case, linear_float, linear_index)
                end
            end
        end
    end
end

@testset "Float32 mini reservoir steps" begin
    for physics in (:single_phase, :immiscible_2ph, :bo_spe1,
            :compositional_2ph_3c, :geothermal)
        @testset "$physics" begin
            case = JutulDarcy.setup_mini_wellcase(Val(physics);
                nstep = 1, total_time = f32_step,
                backend = :csr, block_backend = true)
            precond = :cpr
            if physics == :single_phase
                precond = :ka_spai0
            end
            for (linear_float, linear_index) in
                    ((Float32, Int32), (Float64, Int64))
                @testset "linear $linear_float/$linear_index" begin
                    test_float32_reservoir_case(
                        case, linear_float, linear_index;
                        precond = precond)
                end
            end
        end
    end
end
