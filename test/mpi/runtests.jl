using Test, Jutul, JutulDarcy, MPI, PartitionedArrays, HYPRE

MPI.Init()
@assert MPI.Comm_size(MPI.COMM_WORLD) == 2 "Run this test with test/mpi/run.jl (two MPI ranks)."

function setup_mpi_case()
    bar, day = si_units(:bar, :day)
    domain = reservoir_domain(CartesianMesh((16, 1, 1), (160.0, 10.0, 10.0)),
        permeability = 1e-12, porosity = 0.2)
    rho = [1000.0, 100.0]
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = rho)
    injector = setup_well(domain, [1], name = :Injector)
    producer = setup_well(domain, [16], name = :Producer)
    model = setup_reservoir_model(domain, system;
        wells = [injector, producer], split_wells = true, backend = :csr)
    densities = ConstantCompressibilityDensities(
        p_ref = 100bar, density_ref = rho, compressibility = [1e-6/bar, 1e-4/bar])
    replace_variables!(model, PhaseMassDensities = densities)
    parameters = setup_parameters(model)
    state0 = setup_reservoir_state(model, Pressure = 100bar, Saturations = [0.8, 0.2])
    control = Dict(
        :Injector => InjectorControl(TotalRateTarget(0.001), [0.0, 1.0], density = rho[2]),
        :Producer => ProducerControl(BottomHolePressureTarget(90bar))
    )
    forces = setup_reservoir_forces(model; control = control)
    return JutulCase(model, [day, day], forces; state0 = state0, parameters = parameters)
end

@testset "MPI CPR owned CSR rows" begin
    case = setup_mpi_case()
    options = (info_level = -1, timesteps = :none, failure_cuts_timestep = false,
        error_on_incomplete = true, precond = :cpr, rtol = 1e-8,
        linear_solver_arg = (amg_type = :hypre,))
    sim, config = setup_reservoir_simulator(case; mode = :mpi,
        output_path = only(ARGS), options...)

    # Ensure this really exercises a distributed CSR matrix with ghost rows.
    local_sim = PartitionedArrays.getany(sim.storage.simulators)
    system = Jutul.get_simulator_storage(local_sim).LinearizedSystem
    A = JutulDarcy.reservoir_jacobian(system)
    @test sim.backend isa Jutul.MPI_PArrayBackend
    @test A isa Jutul.StaticSparsityMatrixCSR
    @test 0 < sim.storage.nc_process < size(A, 1)
    @test MPI.Allreduce(sim.storage.nc_process, +, MPI.COMM_WORLD) == 16
    @test config[:linear_solver].preconditioner isa CPRPreconditioner{<:Jutul.BoomerAMGPreconditioner}

    raw_result = simulate!(sim, case.dt; forces = case.forces, config = config)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        # Distributed output is consolidated on rank zero. Compare all report
        # steps, checking the counts first so an empty result cannot pass.
        result = JutulDarcy.ReservoirSimResult(case, raw_result)
        reference = simulate_reservoir(setup_mpi_case(); options...,
            tol_cnv = 1e-7, tol_mb = 1e-10)
        @test length(result.states) == length(reference.states) == length(case.dt)
        @test length(result.time) == length(case.dt)
        for (state, ref) in zip(result.states, reference.states)
            @test state[:Pressure] ≈ ref[:Pressure] rtol = 1e-5
            saturation_error = maximum(abs, state[:Saturations] - ref[:Saturations])
            @test saturation_error < 1e-4
        end
        for well in (:Injector, :Producer), quantity in (:bhp, :rate)
            @test result.wells[well, quantity] ≈ reference.wells[well, quantity] rtol = 1e-5
        end
    end
end
