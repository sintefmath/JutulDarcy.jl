using Jutul, JutulDarcy
using JLArrays
using SparseArrays
using Test

@testset "Single-phase reservoir on a KA backend" begin
    grid = CartesianMesh((4, 1), (4.0, 1.0))
    state0, model, parameters, forces, timesteps = get_test_setup(
        grid,
        case_name = "single_phase_simple",
        context = ParallelCSRContext(1),
        timesteps = [0.1]
    )

    reference_simulator = Simulator(model, state0 = state0, parameters = parameters)
    reference, = simulate!(reference_simulator, timesteps;
        forces = forces, info_level = -1)

    cpu_simulator = Simulator(model, state0 = state0, parameters = parameters)
    simulator = transfer_to_backend(cpu_simulator, JLBackend())
    @test simulator.storage.state.Pressure isa JLArray
    @test simulator.storage.primary_variables.Pressure === simulator.storage.state.Pressure
    @test simulator.storage.LinearizedSystem.jac_buffer ===
        nonzeros(simulator.storage.LinearizedSystem.jac)
    @test cpu_simulator.storage.state.Pressure isa Vector

    states, = simulate!(simulator, timesteps; forces = forces, info_level = -1)
    @test Array(states[end][:Pressure]) ≈ reference[end][:Pressure] rtol = 1e-10
end

@testset "Two-phase reservoir on a KA backend" begin
    grid = CartesianMesh((4, 1), (4.0, 1.0))
    state0, model, parameters, forces, timesteps = get_test_setup(
        grid,
        case_name = "two_phase_simple",
        context = ParallelCSRContext(1),
        timesteps = [0.1]
    )

    reference_simulator = Simulator(model, state0 = state0, parameters = parameters)
    reference, = simulate!(reference_simulator, timesteps;
        forces = forces, info_level = -1)

    cpu_simulator = Simulator(model, state0 = state0, parameters = parameters)
    simulator = transfer_to_backend(cpu_simulator, JLBackend())
    states, = simulate!(simulator, timesteps; forces = forces, info_level = -1)

    @test Array(states[end][:Pressure]) ≈ reference[end][:Pressure] rtol = 1e-10
    @test Array(states[end][:Saturations]) ≈ reference[end][:Saturations] rtol = 1e-10
end
