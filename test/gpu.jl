using Jutul, JutulDarcy, CUDA, SparseArrays, Test
import Adapt

struct PreconvertedLaunchProbe{A}
    output::A
end

function (probe::PreconvertedLaunchProbe)(index, offset)
    probe.output[index] = index + offset
    return nothing
end

function Adapt.adapt_structure(
        ::CUDA.KernelAdaptor, ::PreconvertedLaunchProbe)
    error("Preconverted launch attempted to adapt its callable")
end

function indirection_map_cuda_kernel!(out, map)
    index = Int(CUDA.threadIdx().x)
    if index <= length(map)
        total = zero(eltype(out))
        for value in map[index]
            total += value
        end
        out[index] = total
    end
    return nothing
end

function solve_bl_lsolve(; nx = 10, ny = 1, nstep = nx*ny, lsolve = missing, backend = :csr, step_limit = nothing, kwarg...)
    time = 1.0
    T = time
    tstep = repeat([T/nstep], nstep)
    mesh = CartesianMesh((nx, ny))
    domain = reservoir_domain(mesh)
    nc = number_of_cells(domain)
    timesteps = tstep*3600*24
    bar = 1e5
    p0 = 100*bar
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = [100.0, 100.0])
    model = setup_reservoir_model(domain, sys, extra_out = false, backend = backend)
    kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0], [0.2, 0.2])
    rmodel = reservoir_model(model)
    c = [1e-3, 1e-3]/bar
    density = ConstantCompressibilityDensities(
        p_ref = 100*bar,
        density_ref = [100.0, 100.0],
        compressibility = c
    )
    replace_variables!(rmodel, RelativePermeabilities = kr, PhaseMassDensities = density)
    tot_time = sum(timesteps)
    pv = pore_volume(domain)
    irate = 500*sum(pv)/tot_time
    src  = SourceTerm(nc, irate, fractional_flow = [0.75, 0.25])
    bc = FlowBoundaryCondition(1, p0/2)
    forces = setup_reservoir_forces(model, sources = src, bc = bc)
    parameters = setup_parameters(model)
    state0 = setup_reservoir_state(model, Pressure = p0, Saturations = [0.25, 0.75])
    if ismissing(lsolve)
        lsolve = select_reservoir_linear_solver(model)
    end
    if !isnothing(step_limit)
        timesteps = timesteps[1:step_limit]
    end
    states, report = simulate(state0, model, timesteps; failure_cuts_timestep = false,
        forces = forces, parameters = parameters, linear_solver = lsolve, kwarg...)
    return states[end][:Reservoir][:Saturations][1, :]
end

if CUDA.functional()
    @testset "IndirectionMap on CUDA" begin
        host_map = Jutul.IndirectionMap(Int[1, 2, 3], Int[1, 3, 4])
        device_map = Adapt.adapt(CUDA.CuArray, host_map)
        @test device_map.vals isa CUDA.CuArray
        @test device_map.pos isa CUDA.CuArray
        output = CUDA.zeros(Int, length(host_map))
        CUDA.@cuda threads=length(host_map) indirection_map_cuda_kernel!(
            output, device_map)
        @test Array(output) == [3, 3]
    end
    @testset "Merged wells with host assembly on CUDA" begin
        for simple_well in (true, false)
            bar = si_unit(:bar)
            grid = CartesianMesh((2, 1, 2), (200.0, 100.0, 40.0))
            domain = reservoir_domain(grid,
                permeability = 0.1*si_unit(:darcy), porosity = 0.2)
            wells = [
                setup_vertical_well(domain, 1, 1,
                    name = :Producer, simple_well = simple_well),
                setup_vertical_well(domain, 2, 1,
                    name = :Injector, simple_well = simple_well),
            ]
            system = ImmiscibleSystem((LiquidPhase(), VaporPhase()),
                reference_densities = [1000.0, 100.0])
            model = setup_reservoir_model(domain, system,
                wells = wells, block_backend = true)
            state0 = setup_reservoir_state(model,
                Pressure = 200bar, Saturations = [0.8, 0.2])
            controls = Dict(
                :Producer => ProducerControl(BottomHolePressureTarget(190bar)),
                :Injector => InjectorControl(BottomHolePressureTarget(210bar),
                    [1.0, 0.0], density = 1000.0),
            )
            if simple_well
                controls[:Producer] = DisabledControl()
            end
            forces = setup_reservoir_forces(model, control = controls)
            if simple_well
                producer = model[:Producer]
                nperf = length(physical_representation(producer).perforations.reservoir)
                forces[:Producer] = setup_forces(producer,
                    mask = PerforationMask(ones(nperf)))
            end
            case = JutulDarcy.merge_similar_wells(JutulCase(model,
                [si_unit(:day)], forces, state0 = state0))
            result = simulate_reservoir(case;
                mode = :ka_cuda, wells_on_device = false, info_level = -1)
            @test length(result.states) == 1
            @test all(isfinite, only(result.states)[:Pressure])
        end
    end
    @testset "Compositional flash on CUDA" begin
        data_path = JutulDarcy.GeoEnergyIO.test_input_file_path(
            "SIMPLE_COMP", "SIMPLE_COMP.DATA")
        case = setup_case_from_data_file(data_path;
            flash_reuse_guess = true,
            flash_stability_bypass = true)[1:1]
        result = simulate_reservoir(case;
            mode = :ka_cuda, info_level = -1)
        @test length(result.states) == 1
        @test all(isfinite, only(result.states)[:Pressure])
    end
    @testset "CUDA backend copies" begin
        backend = CUDA.CUDABackend()
        host = zeros(Float64, 4096)
        device = CUDA.CuArray(collect(Float64, eachindex(host)))
        Jutul.prepare_host_transfer!(backend, host)
        Jutul.backend_copyto!(host, device)
        CUDA.synchronize()
        @test host == collect(Float64, eachindex(host))
        allocated = @allocated begin
            Jutul.backend_copyto!(host, device)
            CUDA.synchronize()
        end
        @test allocated < 4096

        # Hybrid Float32 simulations retain Float64 host well state while
        # their device mirror uses Float32 storage.
        mixed_host = zeros(Float64, 3)
        mixed_device = CUDA.CuArray(Float32[1.25, 2.5, 3.75])
        Jutul.backend_copyto!(mixed_host, mixed_device)
        @test mixed_host == Float64[1.25, 2.5, 3.75]
        Jutul.backend_copyto!(mixed_device, Float64[4.5, 5.25, 6.75])
        @test Array(mixed_device) == Float32[4.5, 5.25, 6.75]
    end
    @testset "EGG Float32 CUDA linearization and solve" begin
        path = JutulDarcy.GeoEnergyIO.test_input_file_path(
            "EGG", "EGG.DATA")
        original = setup_case_from_data_file(path;
            backend = :csr, block_backend = true)[1:1]
        case = JutulCase(original.model,
            [0.01*si_unit(:day)], original.forces[1:1];
            state0 = original.state0,
            parameters = original.parameters)
        simulator, config = setup_reservoir_simulator(case;
            mode = :ka_cuda,
            float_type = Float32, index_type = Int32,
            info_level = -1, timesteps = :none)
        reservoir_model = simulator.model.models[:Reservoir]
        @test typeof(first(JutulDarcy.reference_densities(
            reservoir_model.system))) === Float32
        for phase_table in reservoir_model.secondary_variables.PhaseViscosities.pvt
            @test typeof(first(phase_table.tab).p_ref) === Float32
        end
        @test typeof(first(
            reservoir_model.secondary_variables.RelativePermeabilities.krw
        ).connate) === Float32

        forces = Jutul.preprocess_forces(
            simulator, only(case.forces)).forces
        dt = only(case.dt)
        Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
        Jutul.update_state_dependents!(
            simulator.storage, simulator.model, dt, forces; time = dt)
        Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
        host = simulator.storage.host_evaluation
        pressure_drop_wells = filter(host.keys) do key
            haskey(host.storage[key].state, :ConnectionPressureDrop)
        end
        @test !isempty(pressure_drop_wells)
        for key in pressure_drop_wells
            host_dp = host.storage[key].state.ConnectionPressureDrop
            device_dp = simulator.storage[key].state.ConnectionPressureDrop
            @test eltype(host_dp) == Float64
            @test eltype(device_dp) == Float32
            @test host_dp == Array(device_dp)
        end
        Jutul.update_linearized_system!(
            simulator.storage, simulator.model)

        system = simulator.storage.LinearizedSystem
        @test eltype(system.r_buffer) == Float32
        @test all(isfinite, Array(system.r_buffer))
        for block in system.subsystems
            @test all(matrix -> all(isfinite, matrix),
                Array(nonzeros(block.jac)))
        end
        result = Jutul.linear_solve!(system, config[:linear_solver],
            simulator.model, simulator.storage, dt)
        @test result.iterations > 0
        @test all(isfinite, Array(system.dx_buffer))
    end
    @testset "Preconverted CUDA launch" begin
        output = CUDA.zeros(Int, 4)
        probe = PreconvertedLaunchProbe(CUDA.cudaconvert(output))
        context = KernelAbstractionsContext(CUDA.CUDABackend())
        Jutul.KernelExecution.launch_preconverted_threaded_loop(
            probe, length(output), context, 10)
        Jutul.synchronize(context)
        @test Array(output) == collect(11:14)
    end
    do_solve(x) = solve_bl_lsolve(
        lsolve = x,
        info_level = -1
    )
    @testset "KA heterogeneous deck PVT" begin
        grid = CartesianMesh((4, 1), (4.0, 1.0))
        state0, model, parameters, forces, timesteps = get_test_setup(
            grid,
            case_name = "two_phase_simple",
            context = ParallelCSRContext(1),
            timesteps = [0.1]
        )
        water = JutulDarcy.PVTW(ConstMuBTable(
            1.0e7, 1.0, 1.0e-9, 1.0e-3, 0.0))
        oil = JutulDarcy.PVDO(MuBTable(
            [1.0e7, 2.0e7], [1.0, 1.1], [2.0e-3, 2.2e-3]))
        pvt = (water, oil)
        replace_variables!(model,
            PhaseMassDensities = DeckPhaseMassDensities(pvt),
            PhaseViscosities = DeckPhaseViscosities(pvt))

        cpu_simulator = Simulator(
            model, state0 = state0, parameters = parameters)
        simulator = transfer_to_backend(
            cpu_simulator, CUDA.CUDABackend())
        @test haskey(simulator.storage, :evaluation_state)
        @test haskey(simulator.storage, :evaluation_state0)
        @test isbitstype(typeof(evaluation_state(simulator.storage)))
        @test isbitstype(typeof(evaluation_state0(simulator.storage)))
        @test evaluation_state(simulator.storage).Pressure isa
            CUDA.CuDeviceArray
        @test Adapt.adapt(
            CUDA.KernelAdaptor(), evaluation_state(simulator.storage)) ===
            evaluation_state(simulator.storage)
        cached_simulator = transfer_to_backend(
            Simulator(model, state0 = state0, parameters = parameters),
            CUDA.CUDABackend(); reduce_memory = false)

        function assemble(simulator)
            local_forces = Jutul.preprocess_forces(
                simulator, deepcopy(forces)).forces
            dt = only(timesteps)
            Jutul.update_before_step!(
                simulator, dt, local_forces; time = 0.0)
            Jutul.update_state_dependents!(simulator.storage,
                simulator.model, dt, local_forces; time = dt)
            Jutul.update_linearized_system!(
                simulator.storage, simulator.model)
            system = simulator.storage.LinearizedSystem
            return Array(system.r_buffer), Array(nonzeros(system.jac))
        end

        fused_residual, fused_jacobian = assemble(simulator)
        cached_residual, cached_jacobian = assemble(cached_simulator)
        @test all(isfinite, fused_residual)
        @test fused_residual == cached_residual
        @test fused_jacobian == cached_jacobian
    end
    @testset "SPE1 hybrid KA CUDA assembly" begin
        spe1 = JutulDarcy.GeoEnergyIO.test_input_file_path(
            "SPE1", "SPE1.DATA")
        case = setup_case_from_data_file(spe1;
            block_backend = false)[1:1]
        for well_name in (:PROD, :INJ)
            fill!(case.parameters[well_name][:PerforationGravityDifference],
                1.0)
        end
        simulator, = setup_reservoir_simulator(case;
            mode = :ka_cuda,
            info_level = -1,
            linear_solver = nothing,
            timesteps = :none)

        device_cross_terms = filter(simulator.storage.cross_terms) do storage
            storage isa Jutul.KernelExecution.PreparedCrossTermStorage
        end
        @test !isempty(device_cross_terms)
        @test all(device_cross_terms) do storage
            isbitstype(typeof(getfield(storage, :evaluation)))
        end

        @test Jutul.group_execution_mode(simulator.model, :Reservoir) ==
            SolveFullyOnDevice
        @test Jutul.group_execution_mode(simulator.model, :Facility) ==
            AssembleOnDevice
        @test simulator.storage.host_evaluation.keys ==
            (:PROD, :INJ, :Facility)

        if case.forces isa AbstractVector
            forces = only(case.forces)
        else
            forces = case.forces
        end
        forces = deepcopy(forces)
        well = physical_representation(case.model.models.PROD.domain)
        mask = PerforationMask(ones(length(well.perforations.reservoir)))
        forces[:PROD] = setup_forces(case.model.models.PROD, mask = mask)
        forces = Jutul.preprocess_forces(simulator, forces).forces
        @test forces[:PROD].mask.values isa CUDA.CuArray
        dt = only(case.dt)
        host = simulator.storage.host_evaluation
        Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
        Jutul.update_state_dependents!(
            simulator.storage, simulator.model, dt, forces; time = dt)
        Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
        for well_name in (:PROD, :INJ)
            host_dp = host.storage[well_name].state.ConnectionPressureDrop
            backend_dp = Array(
                simulator.storage[well_name].state.ConnectionPressureDrop)
            @test any(value -> !iszero(value), host_dp)
            @test backend_dp == host_dp
        end
        Jutul.update_state_dependents!(
            simulator.storage, simulator.model, dt, forces; time = dt)
        Jutul.update_linearized_system!(simulator.storage, simulator.model)

        system = simulator.storage.LinearizedSystem
        @test system.r_buffer isa CUDA.CuArray
        @test all(isfinite, Array(system.r_buffer))
        if system isa Jutul.MultiLinearizedSystem
            @test all(block -> all(isfinite, Array(nonzeros(block.jac))),
                system.subsystems)
        else
            @test all(isfinite, Array(nonzeros(system.jac)))
        end
    end
    @testset "SPE1 KA CUDA adjoint solve" begin
        spe1 = JutulDarcy.GeoEnergyIO.test_input_file_path(
            "SPE1", "SPE1.DATA")
        case = setup_case_from_data_file(spe1;
            block_backend = true)[1:1]
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
                mode = :ka_cuda,
                linear_solver = nothing,
                timesteps = :none,
            ))
        value, gradient = problem()
        simulator = problem.cache[:simulator]
        storage = problem.cache[:storage]

        @test isfinite(value)
        @test all(isfinite, gradient)
        @test simulator.storage.host_evaluation.keys ==
            (:PROD, :INJ, :Facility)
        @test storage.forward.storage.LinearizedSystem.r_buffer isa CUDA.CuArray
        @test storage.backward.storage.LinearizedSystem.r_buffer isa CUDA.CuArray
        @test storage.parameter.storage.LinearizedSystem.r_buffer isa Vector
        @test storage.state0_buf isa CUDA.CuArray
        @test storage.dstate0 isa Vector
    end
    @testset "SimulationModel" begin
        krylov_cpu = GenericKrylov(:bicgstab, preconditioner = ILUZeroPreconditioner())
        s_cpu = do_solve(krylov_cpu)

        for T in [Float32, Float64]
            for solver in [:bicgstab, :gmres]
                krylov_cu = JutulDarcy.CUDAReservoirKrylov(solver, Float_t = T)
                s_cu = do_solve(krylov_cu)
                @test s_cpu ≈ s_cu
            end
        end
    end
    @testset "SPE9 KA CPR and CPRW" begin
        spe9 = JutulDarcy.GeoEnergyIO.test_input_file_path(
            "SPE9", "SPE9.DATA")
        case = setup_case_from_data_file(spe9;
            block_backend = true)[1:1]
        linear_solver_arg = (
            amg_type = :hmis,
            smoother_type = :ka_ilu0,
            update_interval = :once,
            update_interval_partial = :iteration,
            partial_update = true,
            amg_arg = (reuse = :operators,),
            max_iterations = 100
        )
        for precond in (:cpr, :cprw)
            result = simulate_reservoir(case;
                mode = :ka_cuda,
                precond = precond,
                linear_solver_arg = linear_solver_arg,
                failure_cuts_timestep = false,
                info_level = -1)
            @test length(result.states) == 1
            pressure = only(result.states)[:Pressure]
            @test pressure isa Vector
            @test all(isfinite, pressure)
        end
    end
    function spe1_gpu_compare(ref, cusolve)
        @test ref.wells[:PROD][:grat] ≈ cusolve.wells[:PROD][:grat] rtol = 1e-2
        @test ref.wells[:PROD][:bhp] ≈ cusolve.wells[:PROD][:bhp] rtol = 1e-2
        newtons = report_stats(ref.result.reports).newtons
        newtons_cu = report_stats(cusolve.result.reports).newtons
        @test newtons_cu < newtons + 5
        linits = report_stats(ref.result.reports).linear_iterations
        linits_cu = report_stats(cusolve.result.reports).linear_iterations
        @test linits_cu < linits + 20
    end
    @testset "High-level CUDA" begin
        spe1_pth = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
        case = setup_case_from_data_file(spe1_pth)
        sim_kwarg = (info_level = -1, failure_cuts_timestep = false)
        @testset "CuSPARSE-ILU(0)" begin
            res_ilu = simulate_reservoir(case; precond = :ilu0, sim_kwarg...)
            res_cuilu = simulate_reservoir(case; linear_solver_backend = :cuda, precond = :ilu0, sim_kwarg...);
            spe1_gpu_compare(res_ilu, res_cuilu)
        end
        @testset "KA assembly with CuSPARSE-ILU(0)" begin
            ka_case = setup_case_from_data_file(spe1_pth;
                block_backend = true)[1:1]
            result = simulate_reservoir(ka_case;
                mode = :ka_cuda,
                linear_solver_backend = :cuda,
                precond = :ilu0,
                sim_kwarg...)
            @test length(result.states) == 1
            @test all(isfinite, only(result.states)[:Pressure])
        end
        if Sys.islinux()
            @testset "AMGX-CPR" begin
                using AMGX
                res_cpr = simulate_reservoir(case; precond = :cpr, sim_kwarg...)
                res_cucpr = simulate_reservoir(case; linear_solver_backend = :cuda, precond = :cpr, sim_kwarg...);
                spe1_gpu_compare(res_cpr, res_cucpr)
                ka_case = setup_case_from_data_file(spe1_pth;
                    block_backend = true)[1:1]
                ka_amgx = simulate_reservoir(ka_case;
                    mode = :ka_cuda,
                    linear_solver_backend = :cuda,
                    precond = :cpr,
                    sim_kwarg...)
                @test length(ka_amgx.states) == 1
                @test all(isfinite, only(ka_amgx.states)[:Pressure])
            end
        end
    end
end
