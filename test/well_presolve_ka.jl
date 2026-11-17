using Test, Jutul, JutulDarcy, CUDA
using Jutul.KernelExecution: KernelAbstractions

function simulate_ka_well_presolve(case; mode=:default, wells_on_device=false,
        single_precision_linear_solver=false)
    solver_options = mode == :default ? (linear_solver=nothing,) : (
        precond=:cprw,
        linear_solver_arg=(amg_type=:pmis, smoother_type=:ka_dilu),
    )
    sim, cfg = setup_reservoir_simulator(case; mode=mode,
        wells_on_device=wells_on_device,
        single_precision_linear_solver=single_precision_linear_solver,
        presolve_wells=true, well_acceptance_factor=1.0,
        well_outer_iterations_limit=Inf,
        info_level=-1, timesteps=:none, max_timestep_cuts=0,
        max_nonlinear_iterations=25, failure_cuts_timestep=false,
        solver_options...)
    result = simulate!(sim, case.dt; forces=case.forces,
        state0=case.state0, parameters=case.parameters, config=cfg)
    @test length(result.states) == length(case.dt)
    _, presolve = Jutul.get_prepare_step_handler(sim)
    @test any(pair -> !isnothing(pair[2].F), presolve.groups_and_linear_solvers)
    return result.states, sim
end

@testset "Well presolving on KA backends" begin
    case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph;
        block_backend=true, nstep=2, total_time=2*si_unit(:day))
    reference, = simulate_ka_well_presolve(case)
    modes = CUDA.functional() ? (:ka, :ka_cuda) : (:ka,)
    CUDA.functional() && CUDA.allowscalar(false)
    for mode in modes
        precisions = mode == :ka_cuda ? (false, true) : (false,)
        for wells_on_device in (false, true), single_precision in precisions
            @testset "$mode, device wells=$wells_on_device, Float32 solve=$single_precision" begin
                states, sim = simulate_ka_well_presolve(case; mode=mode,
                    wells_on_device=wells_on_device,
                    single_precision_linear_solver=single_precision)
                for (state, expected) in zip(states, reference)
                    @test state[:Reservoir][:Pressure] ≈
                        expected[:Reservoir][:Pressure] rtol=1e-5
                    @test state[:Reservoir][:Saturations] ≈
                        expected[:Reservoir][:Saturations] atol=1e-5
                end
                if mode == :ka_cuda
                    _, presolve = Jutul.get_prepare_step_handler(sim)
                    @test all(presolve.groups_and_linear_solvers) do (group, _)
                        sim.storage.LinearizedSystem[group,group].r_buffer isa CUDA.CuArray
                    end
                end
            end
        end
    end
end
