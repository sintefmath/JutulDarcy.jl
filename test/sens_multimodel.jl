using Jutul, JutulDarcy, Test

function well_test_objective(model, state)
    if true
        q = state[:Facility][:TotalSurfaceMassRate]
        obj = 2.0*q[1] + 0.5*q[2]
    else
        obj = state[:Reservoir][:Pressure][end]
    end
    return obj
end

function well_test_objective_vec(model, state)
    return [well_test_objective(model, state), well_test_objective(model, state)]
end

function solve_adjoint_forward_test_system()
    states, reports, setup = JutulDarcy.simulate_mini_wellcase(:immiscible_2ph, block_backend = false, general_ad = true)
    return (setup[:model], setup[:state0], states, reports, setup[:parameters], setup[:forces])
end

function test_basic_adjoint()
    model, state0, states, reports, param, forces = solve_adjoint_forward_test_system()
    # Scalar mode - we test the gradient of the scalar objective against the numerical version
    # Define objective
    G = (model, state, dt, step_no, forces) -> well_test_objective(model, state)
    grad_adj = solve_out_of_place(model, state0, states, param, reports, G, forces)
    # Check against numerical gradient
    grad_num = Jutul.solve_numerical_sensitivities(model, states, reports, G, 
                    forces = forces, state0 = state0, parameters = param,
                    epsilon = 1e-4)
    for k in keys(grad_num)
        for ki in keys(grad_num[k])
            adj = grad_adj[k][ki]
            num = grad_num[k][ki]
            @info "$k.$ki" adj num
            # @test isapprox(adj, num, atol = 1e-4)
        end 
    end
end

function test_optimization_gradient(; use_scaling = true, use_log = false)
    model, state0, states, reports, param, forces = solve_adjoint_forward_test_system()
    dt = report_timesteps(reports)
    ϵ = 1e-4
    num_tol = 1e-3

    G = (model, state, dt, step_no, forces) -> well_test_objective(model, state)

    active = nothing
    active = Dict(:Reservoir => [:FluidVolume, :Transmissibilities], :Injector => [:FluidVolume, :WellIndices], :Producer => [:FluidVolume, :WellIndices], :Facility => [])
    # active = Dict(:Reservoir => [:FluidVolume, :Transmissibilities], :Injector => [], :Producer => [], :Facility => [])
    # active = Dict(:Reservoir => [:Transmissibilities], :Injector => [], :Producer => [], :Facility => [])
    active = Dict(:Reservoir => [:PhaseViscosities], :Injector => [], :Producer => [], :Facility => [])
    cfg = optimization_config(model, param, active, use_scaling = use_scaling, rel_min = 0.5, rel_max = 2.0)
    if use_log
        for (k, v) in cfg
            for (ki, vi) in v
                vi[:scaler] = :log
            end
        end
    end
    for (k, v) in cfg
        for (ki, vi) in v
            @info "$k.$ki" vi
        end
    end
    tmp = setup_parameter_optimization(model, state0, param, dt, forces, G, cfg, param_obj = true, print = false);
    F_o, dF_o, F_and_dF, x0, lims, data = tmp
    @info "" data[:mapper]
    for (k, v) in data[:mapper]
        for (ki, vi) in v
            @info "$k.$ki" vi
        end
    end
    # Evaluate gradient first to initialize
    F0 = F_o(x0)
    # This interface is only safe if F0 was called with x0 first.
    dF_initial = dF_o(similar(x0), x0)

    dF_num = similar(dF_initial)
    num_grad!(dF_num, x0, ϵ, F_o)
    @info "" hcat(dF_num, dF_initial)

    # Check around initial point
    @test isapprox(dF_num, dF_initial, rtol = num_tol)
    # Perturb the data in a few different directions and verify
    # the gradients there too. Use the F_and_dF interface, that
    # computes gradients together with the objective
    for delta in [1.05, 0.85, 0.325, 1.55]
        x_mod = delta.*x0
        dF_mod = similar(dF_initial)
        F_and_dF(NaN, dF_mod, x_mod)
        num_grad!(dF_num, x_mod, ϵ, F_o)
        @test isapprox(dF_num, dF_mod, rtol = num_tol)
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


# test_basic_adjoint()
# test_optimization_gradient(use_scaling = false)
test_optimization_gradient(use_scaling = true)

error()
##
@testset "simple adjoint sensitivities" begin
    @testset "adjoint" begin
        for scalar_obj in true # [true, false]
            for in_place in [false, true]
                # Test single step since it hits less of the code
                # test_basic_adjoint(dt = [1.0], in_place = in_place, scalar_obj = scalar_obj)
                # Test with multiple time-steps
                test_basic_adjoint(in_place = in_place, scalar_obj = scalar_obj)
            end
        end
    end
    @testset "optimization interface" begin
        for use_log in [true, false]
            for use_scaling in [true, false]
                # test_optimization_gradient(use_scaling = use_scaling, use_log = use_log)
            end
        end
    end    
end

