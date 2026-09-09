using Jutul, JutulDarcy
using JLArrays
using SparseArrays
using Test

@testset "Convergence reductions on a KA backend" begin
    residual = [1.0 -2.0 3.0; -4.0 5.0 -6.0]
    pore_volume = [2.0, 4.0, 5.0]
    density = [10.0 12.0 14.0; 20.0 22.0 24.0]
    shrinkage = [0.9 0.8 0.7; 0.6 0.5 0.4]
    reference_density = (1000.0, 800.0)
    dt = 0.25
    phases = Val(2)

    cpu_cnv_mb = JutulDarcy.cnv_mb_errors(
        residual, pore_volume, density, dt, phases)
    cpu_bo = JutulDarcy.cnv_mb_errors_bo(
        residual, pore_volume, shrinkage, dt, reference_density, phases)

    context = KernelAbstractionsContext(JLBackend())
    device_cnv_mb = JutulDarcy.cnv_mb_errors(
        JLArray(residual), JLArray(pore_volume), JLArray(density),
        dt, phases, context)
    device_bo = JutulDarcy.cnv_mb_errors_bo(
        JLArray(residual), JLArray(pore_volume), JLArray(shrinkage),
        dt, reference_density, phases, context)

    @test collect(device_cnv_mb[1]) ≈ collect(cpu_cnv_mb[1])
    @test collect(device_cnv_mb[2]) ≈ collect(cpu_cnv_mb[2])
    @test collect(device_bo[1]) ≈ collect(cpu_bo[1])
    @test collect(device_bo[2]) ≈ collect(cpu_bo[2])
end

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
