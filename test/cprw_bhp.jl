using Jutul, JutulDarcy, Test, SparseArrays, CUDA
using Jutul.KernelExecution: KernelAbstractions
using LinearAlgebra: dot

function cprw_bhp_test_case(; merged = false)
    bar = si_unit(:bar)
    grid = CartesianMesh((3, 1, 2), (300.0, 100.0, 40.0))
    domain = reservoir_domain(grid,
        permeability = 0.1*si_unit(:darcy), porosity = 0.2)
    wells = [setup_vertical_well(domain, i, 1; name, simple_well = true)
        for (i, name) in enumerate((:Producer, :Injector, :Disabled))]
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase());
        reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, system; wells,
        backend = :csr, block_backend = true)
    state0 = setup_reservoir_state(model;
        Pressure = 200bar, Saturations = [0.8, 0.2])
    controls = Dict(
        :Producer => ProducerControl(BottomHolePressureTarget(190bar)),
        :Injector => InjectorControl(TotalRateTarget(1e-4),
            [1.0, 0.0]; density = 1000.0),
        :Disabled => DisabledControl(),
    )
    next_controls = copy(controls)
    next_controls[:Producer] = ProducerControl(TotalRateTarget(-1e-4))
    next_controls[:Injector] = InjectorControl(BottomHolePressureTarget(210bar),
        [1.0, 0.0]; density = 1000.0)
    forces = [setup_reservoir_forces(model; control)
        for control in (controls, next_controls)]
    case = JutulCase(model, fill(0.1*si_unit(:day), 2), forces; state0)
    if merged
        case = JutulDarcy.merge_similar_wells(case)
    end
    return case, controls, next_controls
end

function cprw_bhp_test_matrix(reference, backend)
    host = Jutul.KAPreconditioners.csr_matrix(sparse(reference))
    if isnothing(backend)
        return host
    end
    context = KernelAbstractionsContext(backend; float_type = eltype(reference))
    adapt(x) = JutulDarcy.Adapt.adapt(context, x)
    return Jutul.StaticSparsityMatrixCSR(
        adapt(nonzeros(host)), adapt(Jutul.colvals(host)), adapt(host.rowptr),
        size(host)..., backend)
end

@testset "CPRW BHP rows on CPU and device backends" begin
    backends = Any[nothing, KernelAbstractions.CPU()]
    if CUDA.functional()
        CUDA.allowscalar(false)
        push!(backends, CUDA.CUDABackend())
    end
    for merged in (false, true), T in (Float32, Float64), backend in backends
        @testset "merged=$merged, $T, $backend" begin
            case, controls, next_controls = cprw_bhp_test_case(; merged)
            # Six reservoir rows, followed by three independently controlled wells.
            map = if merged
                (wells = [:SimpleWells], well_ranges = [1:3], np = 9)
            else
                (wells = [:Producer, :Injector, :Disabled],
                    well_ranges = [1:1, 2:2, 3:3], np = 9)
            end
            reference = reshape(T.(1:81), 9, 9)
            cfg = JutulDarcy.WellGroupConfiguration(keys(controls), controls)
            storage = Dict(:Facility => (state = (WellGroupConfiguration = cfg,),))
            for (operating, bhp_row) in ((controls, 7), (next_controls, 8))
                # Requested controls stay unchanged when an operating limit switches BHP.
                merge!(cfg.operating_controls, operating)
                matrix = cprw_bhp_test_matrix(reference, backend)
                @test eltype(matrix) === T
                JutulDarcy.cprw_fix_bhp_controlled_wells!(
                    matrix, map, case.model, storage)
                expected = copy(reference)
                expected[bhp_row, :] .= 0
                expected[bhp_row, bhp_row] = reference[bhp_row, bhp_row]
                @test Matrix(JutulDarcy.cpr_host_matrix(matrix)) == expected
            end
            # Missing controls and models without a facility leave the matrix intact.
            empty!(cfg.operating_controls)
            for state in (storage, Dict())
                matrix = cprw_bhp_test_matrix(reference, backend)
                JutulDarcy.cprw_fix_bhp_controlled_wells!(matrix, map, case.model, state)
                @test Matrix(JutulDarcy.cpr_host_matrix(matrix)) == reference
            end
        end
    end
end

@testset "CPR and CPRW adjoint pressure stages" begin
    case = JutulDarcy.setup_mini_wellcase(:compositional_2ph_3c; nstep = 1, backend = :csr)
    (; model, state0, parameters, forces) = case
    # Different weights in each cell expose row/column indexing errors.
    state0[:Reservoir][:Pressure] .*= [0.9, 1.05, 1.2]
    state0[:Reservoir][:OverallMoleFractions][:, 3] .= [0.55, 0.30, 0.15]
    dt = only(case.dt)
    for use_ka in (false, true)
        @testset "KA=$use_ka" begin
            forward = Simulator(model; state0, parameters)
            adjoint = Simulator(Jutul.adjoint_model_copy(model); state0, parameters)
            if use_ka
                forward = transfer_to_backend(forward, KernelAbstractions.CPU())
                adjoint = Jutul.transfer_adjoint_simulator(adjoint, forward.model)
            end
            for sim in (forward, adjoint)
                Jutul.adjoint_reassemble!(sim, state0, state0, dt, forces, dt)
            end
            n = length(Jutul.vector_residual(forward.storage.LinearizedSystem))
            x = collect(range(0.1, 1.0; length = n))
            y = cos.(1:n)
            for variant in (:cpr, :cprw)
                pressures = Matrix{Float64}[]
                corrections = Vector{Float64}[]
                for (mode, sim, residual) in ((:forward, forward, x), (:adjoint, adjoint, y))
                    cpr = select_reservoir_linear_solver(sim.model, variant; mode).preconditioner
                    JutulDarcy.update_cpr_internals!(cpr, sim.storage.LinearizedSystem,
                        sim.model, sim.storage, Jutul.progress_recorder(sim),
                        Jutul.default_executor(), Float64)
                    s = cpr.storage
                    pressure = Matrix(s.A_p)
                    push!(pressures, pressure)
                    JutulDarcy.update_p_rhs!(s.r_p, residual, s.number_of_components,
                        s.block_size, s.w_p, s.A_p, mode)
                    dp = pressure \ s.r_p
                    correction = zeros(n)
                    if mode == :forward
                        JutulDarcy.increment_pressure!(correction, dp, s.block_size)
                    else
                        JutulDarcy.increment_adjoint_pressure!(correction, dp, s.block_size, s.w_p)
                    end
                    push!(corrections, correction)
                end
                @test pressures[2] ≈ pressures[1]'
                @test dot(y, corrections[1]) ≈ dot(x, corrections[2])
            end
            @test Jutul.select_linear_solver(forward.model; mode = :adjoint).preconditioner.variant == :cprw
        end
    end
end

@testset "CPRW simulations with changing BHP controls" begin
    options = (info_level = -1, timesteps = :none,
        error_on_incomplete = true, max_timestep_cuts = 0,
        nonlinear_tolerance = 1e-7, tol_cnv = 1e-7,
        tol_mb = 1e-9, tol_cnv_well = 1e-7)
    policies = [
        (mode = :default,),
        (mode = :ka, ka_backend = KernelAbstractions.CPU()),
    ]
    if CUDA.functional()
        CUDA.allowscalar(false)
        push!(policies, (mode = :ka_cuda,))
    end
    for merged in (false, true)
        case, = cprw_bhp_test_case(; merged)
        reference = simulate_reservoir(case; precond = :ilu0, options...)
        for policy in policies
            @testset "merged=$merged, $policy" begin
                result = simulate_reservoir(case; precond = :cprw, policy..., options...)
                @test length(result.states) == length(case.dt)
                for (actual, expected) in zip(result.states, reference.states)
                    @test actual[:Pressure] ≈ expected[:Pressure] rtol = 1e-7
                    @test actual[:Saturations] ≈ expected[:Saturations] rtol = 1e-7
                end
                cpr = result.extra[:config][:linear_solver].preconditioner
                pressure = Matrix(JutulDarcy.cpr_host_matrix(cpr.storage.A_p))
                @test !iszero(pressure[8, 8])
                @test all(iszero, pressure[8, [1:7; 9]])
                @test any(!iszero, pressure[7, 1:6])
                @test result.wells[:Producer][:bhp][1] ≈ 190*si_unit(:bar)
                @test result.wells[:Injector][:bhp][2] ≈ 210*si_unit(:bar)
            end
        end
    end
end

@testset "CPRW linear solve with CUDA wells" begin
    if CUDA.functional()
        for merged in (false, true)
            case, = cprw_bhp_test_case(; merged)
            simulator, config = setup_reservoir_simulator(case;
                mode = :ka_cuda, wells_on_device = true, precond = :cprw,
                info_level = -1, timesteps = :none)
            forces = Jutul.preprocess_forces(simulator, case.forces[end]).forces
            dt = last(case.dt)
            Jutul.update_before_step!(simulator, dt, forces; time = 0.0)
            Jutul.update_state_dependents!(
                simulator.storage, simulator.model, dt, forces; time = dt)
            Jutul.update_linearized_system!(simulator.storage, simulator.model)
            system = simulator.storage.LinearizedSystem
            @test system.r_buffer isa CUDA.CuArray
            result = Jutul.linear_solve!(system, config[:linear_solver],
                simulator.model, simulator.storage, dt)
            @test result.ok && result.iterations > 0
            @test all(isfinite, Array(system.dx_buffer))
            cpr = config[:linear_solver].preconditioner
            pressure = Matrix(JutulDarcy.cpr_host_matrix(cpr.storage.A_p))
            @test !iszero(pressure[8, 8])
            @test all(iszero, pressure[8, [1:7; 9]])
            @test any(!iszero, pressure[7, 1:6])
        end
    else
        @test_skip CUDA.functional()
    end
end
