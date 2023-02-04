using Jutul, JutulDarcy, Test, LinearAlgebra

function well_test_objective(model, state)
    q = state[:Facility][:TotalSurfaceMassRate]
    return 2.0*q[1] + 0.5*q[2]
end

function solve_adjoint_forward_test_system(; block_backend = true, kwarg...)
    states, reports, setup = JutulDarcy.simulate_mini_wellcase(:immiscible_2ph; kwarg..., block_backend = block_backend, general_ad = false)
    return (setup[:model], setup[:state0], states, reports, setup[:parameters], setup[:forces])
end

function test_optimization_gradient(; use_scaling = true, use_log = false, kwarg...)
    model, state0, states, reports, param, forces = solve_adjoint_forward_test_system(; kwarg...)
    dt = report_timesteps(reports)
    ϵ = 1e-3
    num_tol = 1e-2

    G = (model, state, dt, step_no, forces) -> well_test_objective(model, state)

    active = nothing
    cfg = optimization_config(model, param, active, use_scaling = use_scaling, rel_min = 0.5, rel_max = 2.0)
    if use_log
        for (k, v) in cfg
            for (ki, vi) in v
                vi[:scaler] = :log
            end
        end
    end
    tmp = setup_parameter_optimization(model, state0, param, dt, forces, G, cfg, param_obj = true, print = false);
    F_o, dF_o, F_and_dF, x0, lims, data = tmp
    # Evaluate gradient first to initialize
    F0 = F_o(x0)
    # This interface is only safe if F0 was called with x0 first.
    dF_initial = dF_o(similar(x0), x0)
    dF_num = similar(dF_initial)
    num_grad!(dF_num, x0, ϵ, F_o)

    # Check around initial point
    @test isapprox(dF_num, dF_initial, rtol = num_tol)
    # Perturb the data in a few different directions and verify
    # the gradients there too. Use the F_and_dF interface, that
    # computes gradients together with the objective
    for delta in [1.05, 0.85, 0.325, 1.15]
        x_mod = delta.*x0
        dF_mod = similar(dF_initial)
        F_and_dF(NaN, dF_mod, x_mod)
        num_grad!(dF_num, x_mod, ϵ, F_o)
        err = norm(dF_num - dF_mod)/norm(dF_num)
        @test err <= num_tol
    end
end

function num_grad!(dF_num, x0, ϵ, F_o)
    F_initial = F_o(x0)
    for i in eachindex(dF_num)
        x = copy(x0)
        x[i] += ϵ
        F_perturbed = F_o(x)
        dF_num[i] = (F_perturbed - F_initial)/ϵ
    end
end

function solve_in_place!(model, state0, states, param, dt, G, forces; kwarg...)
    storage = setup_adjoint_storage(model; state0 = state0, parameters = param, kwarg...)
    grad_adj = zeros(storage.n)
    grad_adj = solve_adjoint_sensitivities!(grad_adj, storage, states, state0, dt, G, forces = forces)
    return grad_adj
end

function solve_out_of_place(model, state0, states, param, reports, G, forces; kwarg...)
    grad_adj = solve_adjoint_sensitivities(model, states, reports, G, 
    forces = forces, state0 = state0, parameters = param; kwarg...)
    return grad_adj
end
##
@testset "optimization interface with wells" begin
    for block in [true, false]
        if block
            b = "block"
        else
            b = "scalar"
        end
        @testset "$b" begin
            for ad in [true, false]
                @testset "scaled (linear)" begin
                    test_optimization_gradient(use_scaling = true, general_ad = ad, block_backend = block)
                end
                @testset "scaled (log)" begin
                    test_optimization_gradient(use_scaling = true, use_log = true,  general_ad = ad, block_backend = block)
                end
            end
        end
    end
end
