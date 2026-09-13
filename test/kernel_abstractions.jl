using Jutul, JutulDarcy
using JLArrays
using SparseArrays
using Test

function setup_spe1_ka_case(; block_backend = false)
    spe1 = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
    return setup_case_from_data_file(spe1;
        block_backend = block_backend)[1:1]
end

function setup_spe1_optimization_test(backend)
    case = setup_spe1_ka_case(block_backend = true)
    transmissibilities = copy(
        case.parameters[:Reservoir][:Transmissibilities])
    function setup_optimization_case(parameters, step_info = missing)
        next_case = deepcopy(case)
        next_case.parameters[:Reservoir][:Transmissibilities] =
            transmissibilities .* only(parameters["multiplier"])
        return next_case
    end
    function objective(model, state, dt, step_info, forces)
        pressure = state[:Reservoir][:Pressure]
        value = zero(eltype(pressure))
        for p in pressure
            value += (p/1.0e7)^2
        end
        return value
    end
    parameters = Dict("multiplier" => [1.1])
    dopt = setup_reservoir_dict_optimization(
        parameters, setup_optimization_case; verbose = false)
    free_optimization_parameter!(dopt, "multiplier";
        abs_min = 0.5, abs_max = 2.0)
    problem = JutulDarcy.reservoir_optimization_problem(dopt, objective;
        simulator_arg = (
            mode = :ka,
            ka_backend = backend,
            linear_solver = nothing,
            timesteps = :none,
        ))
    return problem
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
        dt, phases)
    device_bo = JutulDarcy.cnv_mb_errors_bo(
        JLArray(residual), JLArray(pore_volume), JLArray(shrinkage),
        dt, reference_density, phases)

    @test collect(device_cnv_mb[1]) ≈ collect(cpu_cnv_mb[1])
    @test collect(device_cnv_mb[2]) ≈ collect(cpu_cnv_mb[2])
    @test collect(device_bo[1]) ≈ collect(cpu_bo[1])
    @test collect(device_bo[2]) ≈ collect(cpu_bo[2])
end

@testset "KA reservoir solver configuration" begin
    grid = CartesianMesh((2, 1), (2.0, 1.0))
    state0, model, parameters, _, _ = get_test_setup(
        grid,
        case_name = "two_phase_simple",
        context = ParallelCSRContext(matrix_layout = BlockMajorLayout()),
        timesteps = [0.1]
    )
    simulator = transfer_to_backend(
        Simulator(model; state0 = state0, parameters = parameters),
        JLBackend())
    for variant in (:cpr, :cprw)
        solver = select_reservoir_linear_solver(simulator.model, variant;
            smoother_arg = (damping = 0.8,),
            cpr_arg = (weight_scaling = :none,))
        preconditioner = solver.preconditioner
        @test preconditioner isa CPRPreconditioner
        @test preconditioner.variant == variant
        @test preconditioner.weight_scaling == :none
        @test preconditioner.pressure_precond isa Jutul.AMGPreconditioner
        @test preconditioner.pressure_precond.reuse == :memory
        @test preconditioner.pressure_precond.reuse_partial == :operators
        @test preconditioner.update_interval == :ministep
        @test preconditioner.update_interval_partial == :iteration
        @test preconditioner.partial_update
        @test preconditioner.pressure_precond.options.smoother isa
            Jutul.KAPreconditioners.SPAI0
        @test preconditioner.system_precond isa Jutul.KASmootherPreconditioner
        @test preconditioner.system_precond.config isa Jutul.KAPreconditioners.DILU
        @test preconditioner.system_precond.config.damping == 0.8
    end

    solver = select_reservoir_linear_solver(simulator.model, :cpr;
        update_type = :none,
        update_type_partial = :sparsity,
        amg_arg = (smoother_type = :spai0,),
        smoother_type = :ilu0)
    preconditioner = solver.preconditioner
    @test preconditioner.pressure_precond.reuse == :none
    @test preconditioner.pressure_precond.reuse_partial == :sparsity
    @test preconditioner.pressure_precond.options.smoother isa
        Jutul.KAPreconditioners.SPAI0
    @test preconditioner.system_precond.config isa Jutul.KAPreconditioners.ILU0
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
end

@testset "KA reservoir floating-point and index types" begin
    cases = (
        two_phase = JutulDarcy.setup_mini_wellcase(
            Val(:immiscible_2ph); nstep = 1,
            total_time = 0.01*si_unit(:day), backend = :csr,
            block_backend = false),
        compositional = JutulDarcy.setup_mini_wellcase(
            Val(:compositional_2ph_3c); nstep = 1,
            total_time = 0.01*si_unit(:day), backend = :csr,
            block_backend = false, fast_flash = true),
    )
    for (name, case) in pairs(cases)
        group_execution = if name == :compositional
            Dict(:default => Jutul.AssembleOnDevice)
        else
            missing
        end
        simulator, config = setup_reservoir_simulator(case;
            mode = :ka,
            ka_backend = Jutul.KernelExecution.KernelAbstractions.CPU(),
            group_execution = group_execution,
            float_type = Float32,
            index_type = Int32,
            linear_solver = nothing,
            failure_cuts_timestep = false,
            timesteps = :none,
            info_level = -1)

        @test Jutul.float_type(simulator.model.context) === Float32
        @test Jutul.index_type(simulator.model.context) === Int32
        @test all(Jutul.float_type(submodel.context) === Float32
            for submodel in values(simulator.model.models))
        @test all(Jutul.index_type(submodel.context) === Int32
            for submodel in values(simulator.model.models))
        system = simulator.storage.LinearizedSystem
        @test eltype(system.r_buffer) === Float32
        @test eltype(system.jac.nzval) === Float32
        @test eltype(system.jac.rowptr) === Int32
        @test eltype(system.jac.colval) === Int32
        pressure_type = eltype(simulator.storage.Reservoir.state.Pressure)
        @test typeof(Jutul.value(zero(pressure_type))) === Float32

        states, = simulate!(simulator, case.dt;
            forces = case.forces, state0 = case.state0, config = config)
        @test all(isfinite, Array(states[end][:Reservoir][:Pressure]))
    end

    unsupported = try
        setup_reservoir_simulator(cases.two_phase;
            mode = :default, float_type = Float32, index_type = Int32)
        nothing
    catch error
        error
    end
    @test unsupported isa ArgumentError
    @test occursin("only supported for KernelAbstractions",
        sprint(showerror, unsupported))
end

@testset "SPE1 hybrid multimodel on a KA backend" begin
    case = setup_spe1_ka_case()
    simulator, = setup_reservoir_simulator(case;
        mode = :ka,
        ka_backend = JLBackend(),
        info_level = -1,
        linear_solver = nothing,
        timesteps = :none)

    @test simulator.storage.host_evaluation.keys == (:PROD, :INJ, :Facility)
    @test isnothing(simulator.model.groups)
    @test length(simulator.model.group_execution) ==
        length(simulator.model.models)
    @test Jutul.group_execution_mode(simulator.model, :Reservoir) ==
        SolveFullyOnDevice
    @test Jutul.group_execution_mode(simulator.model, :PROD) ==
        AssembleOnDevice
    @test Jutul.group_execution_mode(simulator.model, :Facility) ==
        AssembleOnDevice
    host = simulator.storage.host_evaluation
    @test host.model.models.Facility.domain.well_symbols isa
        Vector{Symbol}
    @test simulator.model.models.Facility.domain.well_symbols isa Vector{Symbol}
    @test !(simulator.model.models.Facility.domain.well_symbols isa Tuple)
    @test simulator.storage.PROD.state.Pressure isa JLArray
    @test simulator.storage.INJ.state.Pressure isa JLArray
    @test host.storage.PROD.state.Pressure isa Vector
    @test host.storage.INJ.state.Pressure isa Vector
    @test simulator.storage.Facility.state.WellGroupConfiguration !== nothing
    @test host.storage.Facility.state.WellGroupConfiguration !== nothing
    @test all(cross_term ->
            cross_term.target_impact_map.entries isa JLArray,
        simulator.storage.cross_terms)
    mask = PerforationMask([1.0, 0.0])
    adapted_mask = Jutul.preprocess_forces(
        simulator, (mask = mask,)).forces.mask
    @test adapted_mask.values isa JLArray
    reset_state = deepcopy(case.state0)
    reset_state[:PROD][:Pressure] .+= 1.0
    Jutul.reset_variables!(simulator, reset_state)
    @test value.(host.storage.PROD.state.Pressure) ≈
        reset_state[:PROD][:Pressure]
    @test value.(Array(simulator.storage.PROD.state.Pressure)) ≈
        reset_state[:PROD][:Pressure]
    Jutul.reset_variables!(simulator, case.state0)

    if case.forces isa AbstractVector
        forces = first(case.forces)
    else
        forces = case.forces
    end
    forces = Jutul.preprocess_forces(simulator, forces).forces
    dt = first(case.dt)
    Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
    Jutul.update_state_dependents!(
        simulator.storage, simulator.model, dt, forces; time = dt)
    Jutul.update_linearized_system!(simulator.storage, simulator.model)

    system = simulator.storage.LinearizedSystem
    @test system isa Jutul.LinearizedSystem
    @test system.r_buffer isa JLArray
    @test all(isfinite, Array(system.r_buffer))
    @test all(isfinite, Array(nonzeros(system.jac)))

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

    if case.forces isa AbstractVector
        forces = only(case.forces)
    else
        forces = case.forces
    end
    forces = Jutul.preprocess_forces(simulator, forces).forces
    dt = only(case.dt)
    Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
    Jutul.update_state_dependents!(
        simulator.storage, simulator.model, dt, forces; time = dt)
    Jutul.update_linearized_system!(simulator.storage, simulator.model)

    system = simulator.storage.LinearizedSystem
    @test system.r_buffer isa JLArray
    @test all(isfinite, Array(system.r_buffer))
    @test all(isfinite, Array(nonzeros(system.jac)))
end

@testset "SPE1 adjoint solve on a KA backend" begin
    problem = setup_spe1_optimization_test(JLBackend())
    objective, gradient = problem()
    simulator = problem.cache[:simulator]
    storage = problem.cache[:storage]

    @test isfinite(objective)
    @test all(isfinite, gradient)
    @test simulator.storage.host_evaluation.keys ==
        (:PROD, :INJ, :Facility)
    @test length(simulator.model.group_execution) ==
        length(simulator.model.models)
    @test storage.forward.storage.LinearizedSystem.r_buffer isa JLArray
    @test storage.backward.storage.LinearizedSystem.r_buffer isa JLArray
    @test storage.parameter.storage.LinearizedSystem.r_buffer isa Vector
    @test storage.state0_buf isa JLArray
    @test storage.dstate0 isa Vector
end
