using Jutul, JutulDarcy
using JLArrays
using SparseArrays
using Test

function setup_spe1_ka_case()
    spe1 = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
    return setup_case_from_data_file(spe1;
        block_backend = false)[1:1]
end

function apply_ka_cpr!(simulator; variant = :cpr)
    system = simulator.storage.LinearizedSystem
    Jutul.prepare_linear_solve!(system)
    preconditioner = CPRPreconditioner(
        LUPreconditioner(), LUPreconditioner();
        strategy = :true_impes,
        variant = variant
    )
    Jutul.update_preconditioner!(preconditioner, system,
        simulator.model.context, simulator.model, simulator.storage,
        ProgressRecorder(), simulator.executor)
    pressure_matrix = preconditioner.pressure_precond.host_matrix
    pressure_values = preconditioner.pressure_precond.host_source_values
    system_matrix = preconditioner.system_precond.host_matrix
    system_values = preconditioner.system_precond.host_source_values
    Jutul.update_preconditioner!(preconditioner, system,
        simulator.model.context, simulator.model, simulator.storage,
        ProgressRecorder(), simulator.executor)
    cache_reused =
        preconditioner.pressure_precond.host_matrix === pressure_matrix &&
        preconditioner.pressure_precond.host_source_values === pressure_values &&
        preconditioner.system_precond.host_matrix === system_matrix &&
        preconditioner.system_precond.host_source_values === system_values
    rhs = Jutul.vector_residual(system)
    increment = similar(rhs)
    Jutul.apply!(increment, preconditioner, rhs)
    return preconditioner, increment, cache_reused
end

@testset "Heterogeneous deck PVT adaptation" begin
    water = JutulDarcy.PVTW(ConstMuBTable(
        1.0e7, 1.0, 1.0e-9, 1.0e-3, 0.0))
    oil = JutulDarcy.PVDO(MuBTable(
        [1.0e7, 2.0e7], [1.0, 1.1], [2.0e-3, 2.2e-3]))
    variable = DeckPhaseMassDensities((water, oil))
    context = KernelAbstractionsContext(JLBackend())

    adapted = JutulDarcy.Adapt.adapt(context, variable)

    @test adapted.pvt[2].tab[1].pressure isa JLArray
    @test adapted.pvt[2].tab[1].shrinkage isa JLArray
    @test adapted.pvt[2].tab[1].viscosity isa JLArray

    pc = Jutul.LinearInterpolant([0.0, 1.0], [0.0, 1.0])
    scaled_pc = JutulDarcy.ScaledCapillaryPressure(
        ((pc, pc), (pc, pc)); regions = [1])
    adapted_pc = JutulDarcy.Adapt.adapt(context, scaled_pc)
    @test adapted_pc.pc[1][1].X isa JLArray
    @test adapted_pc.regions isa JLArray
end

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

    reset_state = deepcopy(state0)
    reset_state[:Pressure] .+= 1.0
    Jutul.reset_variables!(simulator, reset_state)
    @test value.(Array(simulator.storage.state.Pressure)) ≈
        reset_state[:Pressure]
    Jutul.reset_variables!(simulator, state0)

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

    cpr, increment, cache_reused = apply_ka_cpr!(simulator)
    @test all(isfinite, Array(increment))
    @test cache_reused
    @test cpr.storage.A_p.nzval isa JLArray
    @test cpr.pressure_precond.host_matrix !== nothing
    @test cpr.system_precond.host_matrix !== nothing
end

@testset "SPE1 hybrid multimodel on a KA backend" begin
    case = setup_spe1_ka_case()
    simulator, = setup_reservoir_simulator(case;
        mode = :ka,
        ka_backend = JLBackend(),
        info_level = -1,
        linear_solver = nothing,
        timesteps = :none)

    @test !Jutul.multi_model_is_specialized(simulator.model)
    @test simulator.storage isa Jutul.JutulStorage{Nothing}
    @test simulator.storage.cross_terms isa AbstractVector
    @test simulator.storage.host_evaluation.keys == [:PROD, :INJ, :Facility]
    @test simulator.model.groups == [1, 2, 2, 2]
    @test Jutul.group_execution_mode(simulator.model, :Reservoir) ==
        SolveFullyOnDevice
    @test Jutul.group_execution_mode(simulator.model, :PROD) ==
        AssembleOnDevice
    @test Jutul.group_execution_mode(simulator.model, :Facility) ==
        AssembleOnDevice
    host = simulator.storage.host_evaluation
    @test host.model.models.Facility.domain.well_symbols isa
        Vector{Symbol}
    @test simulator.model.models.Facility.domain.well_symbols isa Base.OneTo
    @test !(simulator.model.models.Facility.domain.well_symbols isa Tuple)
    @test simulator.storage.PROD.state.Pressure isa JLArray
    @test simulator.storage.INJ.state.Pressure isa JLArray
    @test host.storage.PROD.state.Pressure isa Vector
    @test host.storage.INJ.state.Pressure isa Vector
    @test simulator.storage.Facility.state.WellGroupConfiguration === nothing
    @test host.storage.Facility.state.WellGroupConfiguration !== nothing
    @test all(cross_term ->
            cross_term.target_impact_map.entries isa JLArray,
        simulator.storage.cross_terms)
    masked_cross_terms = filter(
        cross_term -> haskey(cross_term, :force_buffer),
        simulator.storage.cross_terms)
    @test !isempty(masked_cross_terms)
    @test all(cross_term -> cross_term.force_buffer isa JLArray,
        masked_cross_terms)

    reset_state = deepcopy(case.state0)
    reset_state[:PROD][:Pressure] .+= 1.0
    Jutul.reset_variables!(simulator, reset_state)
    @test value.(host.storage.PROD.state.Pressure) ≈
        reset_state[:PROD][:Pressure]
    @test value.(Array(simulator.storage.PROD.state.Pressure)) ≈
        reset_state[:PROD][:Pressure]
    Jutul.reset_variables!(simulator, case.state0)

    forces = case.forces isa AbstractVector ? first(case.forces) : case.forces
    dt = first(case.dt)
    Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
    Jutul.update_state_dependents!(
        simulator.storage, simulator.model, dt, forces; time = dt)
    Jutul.update_linearized_system!(simulator.storage, simulator.model)

    system = simulator.storage.LinearizedSystem
    @test system isa Jutul.MultiLinearizedSystem
    @test system.r_buffer isa JLArray
    @test all(isfinite, Array(system.r_buffer))
    @test all(block -> all(isfinite, Array(nonzeros(block.jac))),
        system.subsystems)

    tolerances = Jutul.set_default_tolerances(simulator.model)
    converged, error, errors = Jutul.check_convergence(
        simulator.storage,
        simulator.model,
        Dict(:tolerances => tolerances);
        dt = dt,
        extra_out = true
    )
    @test converged isa Bool
    @test isfinite(error)
    @test Set(keys(errors)) == Set(keys(simulator.model.models))

    fill!(system.dx_buffer, 0.0)
    report = Jutul.update_primary_variables!(simulator.storage, simulator.model)
    @test Set(keys(report)) == Set(keys(simulator.model.models))
end

@testset "SPE1 reservoir simulator assembly on a KA backend" begin
    case = setup_spe1_ka_case()
    simulator, = setup_reservoir_simulator(case;
        mode = :ka,
        ka_backend = JLBackend(),
        info_level = -1,
        linear_solver = nothing,
        timesteps = :none)

    forces = case.forces isa AbstractVector ? only(case.forces) : case.forces
    dt = only(case.dt)
    Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
    Jutul.update_state_dependents!(
        simulator.storage, simulator.model, dt, forces; time = dt)
    Jutul.update_linearized_system!(simulator.storage, simulator.model)

    system = simulator.storage.LinearizedSystem
    @test system.r_buffer isa JLArray
    @test all(isfinite, Array(system.r_buffer))
    @test all(block -> all(isfinite, Array(nonzeros(block.jac))),
        system.subsystems)
    cprw, increment, cache_reused = apply_ka_cpr!(simulator; variant = :cprw)
    @test all(isfinite, Array(increment))
    @test cache_reused
    @test cprw.storage.pressure_map isa NamedTuple
    @test size(cprw.storage.A_p, 1) > cprw.storage.pressure_map.ncell
    @test cprw.storage.A_p.nzval isa JLArray
    @test cprw.pressure_precond.host_matrix !== nothing
    @test cprw.system_precond.host_matrix !== nothing
end
