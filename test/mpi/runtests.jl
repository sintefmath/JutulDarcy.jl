using Test, Jutul, JutulDarcy, MPI, PartitionedArrays, HYPRE

MPI.Init()
@assert MPI.Comm_size(MPI.COMM_WORLD) == 2 "Run this test with test/mpi/run.jl (two MPI ranks)."

function setup_mpi_case(; clustered = false, multisegment = false, thermal = false)
    bar, day = si_units(:bar, :day)
    domain = reservoir_domain(CartesianMesh((16, 1, 1), (160.0, 10.0, 10.0)),
        permeability = 1e-12, porosity = 0.2)
    rho = [1000.0, 100.0]
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = rho)
    injector = setup_well(domain, multisegment ? [1, 2] : [1];
        name = :Injector, simple_well = !multisegment)
    producer_cells = multisegment ? (clustered ? [3, 4] : [15, 16]) : [clustered ? 2 : 16]
    producer = setup_well(domain, producer_cells;
        name = :Producer, simple_well = !multisegment)
    model = setup_reservoir_model(domain, system;
        wells = [injector, producer], split_wells = true, backend = :csr, thermal)
    if thermal
        push!(reservoir_model(model).output_variables, :Temperature)
        unique!(reservoir_model(model).output_variables)
    end
    densities = ConstantCompressibilityDensities(
        p_ref = 100bar, density_ref = rho, compressibility = [1e-6/bar, 1e-4/bar])
    replace_variables!(model, PhaseMassDensities = densities)
    parameters = setup_parameters(model)
    state_options = thermal ? (Temperature = 350.0,) : NamedTuple()
    state0 = setup_reservoir_state(model; Pressure = 100bar, Saturations = [0.8, 0.2], state_options...)
    control = Dict(
        :Injector => InjectorControl(TotalRateTarget(0.001), [0.0, 1.0], density = rho[2], temperature = 350.0),
        :Producer => ProducerControl(BottomHolePressureTarget(90bar))
    )
    forces = setup_reservoir_forces(model; control = control)
    if clustered
        next_control = Dict(
            :Injector => InjectorControl(BottomHolePressureTarget(110bar), [0.0, 1.0], density = rho[2], temperature = 350.0),
            :Producer => ProducerControl(TotalRateTarget(-0.001)),
        )
        forces = [forces, setup_reservoir_forces(model; control = next_control)]
    end
    return JutulCase(model, [day, day], forces; state0 = state0, parameters = parameters)
end

@testset "MPI $precond owned CSR rows (clustered=$clustered, multisegment=$multisegment)" for (precond, clustered, multisegment) in ((:cpr, false, false), (:cprw, false, false), (:cprw, true, true))
    # The last case also checks bz != ncomp (temperature is a fine unknown).
    thermal = clustered && multisegment
    case = setup_mpi_case(; clustered, multisegment, thermal)
    options = (info_level = -1, timesteps = :none, failure_cuts_timestep = false,
        error_on_incomplete = true, precond = precond, rtol = 1e-8,
        linear_solver_arg = (amg_type = :hypre,))
    output_path = joinpath(only(ARGS), "$(precond)_$(clustered)_$(multisegment)")
    if clustered
        # Deliberately unequal cell counts and two wells on rank zero. Rank one
        # owns no wells and must still shift its rows and remote column indices.
        partition = JutulDarcy.reservoir_partition(case.model, [fill(1, 6); fill(2, 10)])
        sim = Jutul.PArraySimulator(case, partition; backend = Jutul.MPI_PArrayBackend())
        solver = select_reservoir_linear_solver(case.model, precond;
            rtol = options.rtol, options.linear_solver_arg...,
            cpr_arg = (npre = 0, npost = 1))
        config = simulator_config(sim; linear_solver = solver, output_path,
            info_level = -1, failure_cuts_timestep = false, error_on_incomplete = true,
            timestep_selectors = Any[])
        JutulDarcy.set_default_cnv_mb!(config, sim; tol_cnv = 1e-7, tol_mb = 1e-10)
    else
        sim, config = setup_reservoir_simulator(case; mode = :mpi, output_path, options...)
    end

    # Ensure this really exercises a distributed CSR matrix with ghost rows.
    local_sim = PartitionedArrays.getany(sim.storage.simulators)
    system = Jutul.get_simulator_storage(local_sim).LinearizedSystem
    A = JutulDarcy.reservoir_jacobian(system)
    @test size(eltype(A), 1) == 2 + thermal
    if thermal
        local_model = Jutul.get_simulator_model(local_sim)
        push!(reservoir_model(local_model).output_variables, :Temperature)
        unique!(reservoir_model(local_model).output_variables)
    end
    @test sim.backend isa Jutul.MPI_PArrayBackend
    @test A isa Jutul.StaticSparsityMatrixCSR
    @test 0 < sim.storage.nc_process < size(A, 1)
    @test MPI.Allreduce(sim.storage.nc_process, +, MPI.COMM_WORLD) == 16
    @test config[:linear_solver].preconditioner isa CPRPreconditioner{<:Jutul.BoomerAMGPreconditioner}

    raw_result = simulate!(sim, case.dt; forces = case.forces, config = config)
    _, local_preconditioners = config[:linear_solver].preconditioner
    local_cpr = PartitionedArrays.getany(local_preconditioners)
    @test local_cpr.storage.A_p isa HYPRE.HYPREMatrix
    @test local_cpr.amg_rebuild_tolerance > 0
    @test isnothing(local_cpr.amg_rebuild_storage)
    pressure_rows = local_cpr.storage.A_p.iupper - local_cpr.storage.A_p.ilower + 1
    @test MPI.Allreduce(pressure_rows, +, MPI.COMM_WORLD) == 16 + 2*(precond == :cprw)
    if clustered
        @test pressure_rows == (MPI.Comm_rank(MPI.COMM_WORLD) == 0 ? 8 : 10)
    end
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        # Distributed output is consolidated on rank zero. Compare all report
        # steps, checking the counts first so an empty result cannot pass.
        result = JutulDarcy.ReservoirSimResult(case, raw_result)
        reference = simulate_reservoir(setup_mpi_case(; clustered, multisegment, thermal); options...,
            tol_cnv = 1e-7, tol_mb = 1e-10)
        @test length(result.states) == length(reference.states) == length(case.dt)
        @test length(result.time) == length(case.dt)
        for (state, ref) in zip(result.states, reference.states)
            @test state[:Pressure] ≈ ref[:Pressure] rtol = 1e-5
            saturation_error = maximum(abs, state[:Saturations] - ref[:Saturations])
            @test saturation_error < 1e-4
            if thermal
                @test state[:Temperature] ≈ ref[:Temperature] rtol = 1e-5
            end
        end
        for well in (:Injector, :Producer), quantity in (:bhp, :rate)
            @test result.wells[well, quantity] ≈ reference.wells[well, quantity] rtol = 1e-5
        end
    end
end
