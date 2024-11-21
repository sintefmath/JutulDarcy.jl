using Jutul, JutulDarcy, CUDA, Test

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
        lsolve = reservoir_linsolve(model)
    end
    if !isnothing(step_limit)
        timesteps = timesteps[1:step_limit]
    end
    states, report = simulate(state0, model, timesteps; failure_cuts_timestep = false,
        forces = forces, parameters = parameters, linear_solver = lsolve, kwarg...)
    return states[end][:Reservoir][:Saturations][1, :]
end

if CUDA.functional()
    do_solve(x) = solve_bl_lsolve(
        lsolve = x,
        info_level = -1
    )
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
        if Sys.islinux()
            @testset "AMGX-CPR" begin
                using AMGX
                res_cpr = simulate_reservoir(case; precond = :cpr, sim_kwarg...)
                res_cucpr = simulate_reservoir(case; linear_solver_backend = :cuda, precond = :cpr, sim_kwarg...);
                spe1_gpu_compare(res_cpr, res_cucpr)
            end
        end
    end
end
