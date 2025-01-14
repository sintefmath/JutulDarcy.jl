using Jutul, JutulDarcy, Test, LinearAlgebra, HYPRE

function well_test_objective(model, state)
    q = state[:Facility][:TotalSurfaceMassRate]
    return 2.0*q[1] + 0.5*q[2]
end

function solve_adjoint_forward_test_system(casename; block_backend = true, kwarg...)
    Darcy = si_unit(:darcy)
    states, reports, setup = JutulDarcy.simulate_mini_wellcase(
        casename;
        kwarg...,
        total_time = 1.0*si_unit(:day),
        nstep = 5,
        block_backend = block_backend,
        )
    return (setup[:model], setup[:state0], states, reports, setup[:parameters], setup[:forces])
end

function test_optimization_gradient(casename = :immiscible_2ph; use_scaling = true, use_log = false, parameter_subset = missing, kwarg...)
    model, state0, states, reports, param, forces = solve_adjoint_forward_test_system(casename; kwarg...)
    dt = report_timesteps(reports)
    ϵ = 1e-6
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
    if !ismissing(parameter_subset)
        for (k, v) in cfg
            for (ki, vi) in v
                vi[:active] = ki in parameter_subset
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
    for delta in [1.05, 0.85, 0.325, 1.2]
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

@testset "optimization interface with wells" begin
    for casename in [:compositional_2ph_3c, :immiscible_2ph, :bo_spe1]
        @testset "$casename" begin
            for block in [true, false]
                if block
                    b = "block"
                else
                    b = "scalar"
                end
                @testset "$b" begin
                    for ad in [true, false]
                        @testset "scaled (linear)" begin
                            test_optimization_gradient(casename, use_scaling = true, general_ad = ad, block_backend = block)
                        end
                        @testset "scaled (log)" begin
                            test_optimization_gradient(casename, use_scaling = true, use_log = true,  general_ad = ad, block_backend = block)
                        end
                        @testset "subsets" begin
                            test_optimization_gradient(casename, parameter_subset = (:WellIndices, :Transmissibilities, :StaticFluidVolume))
                            test_optimization_gradient(casename, parameter_subset = (:WellIndices, ))
                        end
                    end
                end
            end
        end
    end
end

@testset "Egg sensitivities" begin
    egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG");
    data_pth = joinpath(egg_dir, "EGG.DATA");

    case_f   = setup_case_from_data_file(data_pth)[1:20];
    result_f = simulate_reservoir(case_f, info_level = -1);

    case_c   = JutulDarcy.coarsen_reservoir_case(case_f, (6, 6, 1), method = :ijk);
    result_c = simulate_reservoir(case_c, info_level = -1, tol_cnv = 1e-6);

    states_f = result_f.result.states;
    model_f  = case_f.model;

    states_c = result_c.result.states;
    model_c  = case_c.model;
    state0_c = setup_state(model_c, case_c.state0);
    param_c  = setup_parameters(model_c)
    forces_c = case_c.forces;
    dt       = case_c.dt;

    bhp = JutulDarcy.BottomHolePressureTarget(1.0);
    wells = collect(keys(JutulDarcy.get_model_wells(case_f)));

    day = si_unit(:day);
    wrat_scale = (1/150)*day;
    orat_scale = (1/80) *day;
    grat_scale = (1/1000)*day;
    bhp_scale  = 1/(20*si_unit(:bar));

    w = [];
    matches = [];
    signs = [];
    sys = reservoir_model(model_f).system;
    wrat = SurfaceWaterRateTarget(-1.0);
    orat = SurfaceOilRateTarget(-1.0);
    grat = SurfaceGasRateTarget(-1.0);

    push!(matches, bhp);
    push!(w, bhp_scale);
    push!(signs, 0);

    for phase in JutulDarcy.get_phases(sys)
        if phase == LiquidPhase()
            push!(matches, orat)
            push!(w, orat_scale)
            push!(signs, -1)

        elseif phase == VaporPhase()
            push!(matches, grat)
            push!(w, grat_scale)
            push!(signs, -1)
        else
            @assert phase == AqueousPhase()
            push!(matches, wrat)
            push!(w, wrat_scale)
            push!(signs, -1)
        end
    end
    o_scale = 1.0 / (sum(dt) * length(wells));
    G = (model_c, state_c, dt, step_no, forces) -> well_mismatch(
        matches, wells, model_f, states_f, model_c, state_c, dt, step_no,
        forces, weights = w, scale = o_scale, signs = signs
    );

    @test Jutul.evaluate_objective(G, model_c, states_c, dt, case_c.forces) > 0.0
    @test Jutul.evaluate_objective(G, model_f, states_f, dt, case_f.forces) ≈ 0.0

    function get_sens(model, state0, parameters, tstep, forces, G)
        sim = Simulator(model, state0 = state0, parameters = parameters)
        states, reports = simulate(sim, tstep, forces=forces, info_level = -1)

        grad_adj = solve_adjoint_sensitivities(
            model, states, reports, G,
            forces = forces,
            state0 = state0,
            parameters = parameters,
            raw_output = false
        )

        grad_numeric = Dict{Symbol,Any}()

        is_multi = isa(model, MultiModel)
        if is_multi
            grad_adj = grad_adj[:Reservoir]
        end
        for k in keys(grad_adj)
            if is_multi
                target = (:Reservoir, k)
            else
                target = k
            end
            grad_numeric[k] = Jutul.solve_numerical_sensitivities(model, states, reports, G, target,
                forces = forces,
                state0 = state0,
                parameters = parameters,
                epsilon = 1e-6
            )
        end
        return grad_adj, grad_numeric
    end

    adj, num = get_sens(model_c, state0_c, param_c, dt, forces_c, G);
    for k in [:Transmissibilities, :FluidVolume, :TwoPointGravityDifference]
        numgrad = num[k]
        adjgrad = adj[k]
        if k == :Transmissibilities
            # Natural scaling of transmissibilities is large in strict SI units (m/s)
            scale = 1e12
        else
            scale = 1.0
        end
        @test norm(adjgrad - numgrad)/scale < 1e-6
    end
end
