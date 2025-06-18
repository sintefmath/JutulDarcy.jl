using Jutul, JutulDarcy, Test, LinearAlgebra
function setup_bl_twoforces(;nc = 100, time = 1.0, nstep = 100)
    T = time
    tstep = repeat([T/nstep], nstep)
    G = get_1d_reservoir(nc)
    nc = number_of_cells(G)
    timesteps = tstep*3600*24 # Convert time-steps from days to seconds

    bar = 1e5
    p0 = 100*bar
    # Define system and realize on grid
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    model = SimulationModel(G, sys)
    kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0], [0.2, 0.2])
    replace_variables!(model, RelativePermeabilities = kr)
    tot_time = sum(timesteps)
    pv = pore_volume(G)
    irate = 500*sum(pv)/tot_time
    src  = SourceTerm(1, irate, fractional_flow = [1.0, 0.0])
    bc = FlowBoundaryCondition(nc, p0/2)
    forces = setup_forces(model, sources = src, bc = bc)

    src2 = SourceTerm(1, 0.99*irate, fractional_flow = [1.0, 0.0])
    force2 = setup_forces(model, sources = src2, bc = bc)
    forces = repeat([forces], length(tstep))
    for i in 1:(nstep÷2)
        forces[i] = force2
    end

    parameters = setup_parameters(model, PhaseViscosities = [1e-3, 5e-3]) # 1 and 5 cP
    state0 = setup_state(model, Pressure = p0, Saturations = [0.7, 0.3])
    return (model, state0, parameters, forces, tstep)
end

function test_force_vectorization(forces, tstep, model)
    unique_forces, to_step = Jutul.unique_forces_and_mapping(forces, tstep)
    for (i, uf) in enumerate(unique_forces)
        x, cfg = Jutul.vectorize_forces(uf, model)
        new_force = Jutul.devectorize_forces(uf, model, x, cfg)
        if model isa SimulationModel
            uf = Dict(:Model => uf)
            new_force = Dict(:Model => new_force)
        end
        for (k, v) in uf
            @testset "$k force set $i" begin
                for fkey in keys(v)
                    @testset "$fkey" begin
                        @test haskey(new_force, k)
                        f_new = new_force[k][fkey]
                        f_old = v[fkey]

                        if f_new isa Dict && f_old isa Dict
                            for (k2, v2) in pairs(f_old)
                                @testset "$k2" begin
                                    @test haskey(f_new, k2)
                                    @test isequal(f_new[k2], v2)
                                end
                            end
                        else
                            @test isequal(f_new, f_old)
                        end
                    end
                end
            end
        end
    end
end

function numerical_diff_forces(model, state0, parameters, forces, tstep, G, eps = 1e-6; eachstep = false)
    function perturb(sim, x, i, ϵ, F, cfg, ix)
        x_delta = copy(x)
        x_delta[i] += ϵ
        new_force = Jutul.devectorize_forces(F, model, x_delta, cfg)
        new_forces = deepcopy(forces)
        for j in ix
            new_forces[j] = new_force
        end
        s, r = simulate!(sim, tstep, state0 = state0, forces = new_forces, parameters = parameters, info_level = -1)
        return Jutul.evaluate_objective(G, model, s, tstep, new_forces)
    end
    dx = Vector{Float64}[]
    sim = Simulator(model, state0 = state0)
    s0, = simulate!(sim, tstep, state0 = state0, forces = forces, parameters = parameters, info_level = -1)
    obj0 = Jutul.evaluate_objective(G, model, s0, tstep, forces)
    unique_forces, to_step = Jutul.unique_forces_and_mapping(forces, tstep, eachstep = eachstep)
    for fno in eachindex(unique_forces)
        F = unique_forces[fno]
        x, cfg = Jutul.vectorize_forces(F, model)
        dx_i = Float64[]
        for i in eachindex(x)
            ϵ = max(1e-18, eps*abs(x[i]))
            stepix = to_step[fno]
            obj_plus = perturb(sim, x, i, ϵ, F, cfg, stepix)
            obj_minus = perturb(sim, x, i, -ϵ, F, cfg, stepix)
            push!(dx_i, (obj_plus - obj_minus)/(2*ϵ))
        end
        push!(dx, dx_i)
    end
    return dx
end

@testset "bc and source force gradients" begin
    model, state0, parameters, forces, tstep = setup_bl_twoforces(nc = 10, nstep = 10)
    test_force_vectorization(forces, tstep, model)

    G = (model, state, dt, step_info, forces) -> dt*(sum(state[:Saturations][1, :] .- 0.5))^2
    dx = numerical_diff_forces(model, state0, parameters, forces, tstep, G)
    states, reports = simulate(state0, model, tstep, forces = forces, parameters = parameters, info_level = -1)
    # Check numerical gradients
    dforces, t_to_f, grad_adj = Jutul.solve_adjoint_forces(model, states, reports, G, forces,
                    state0 = state0, parameters = parameters)

    for i in eachindex(dx, grad_adj)
        @test isapprox(dx[i], grad_adj[i], atol = 1e-3, rtol = 1e-3)
    end
    # Check optimization interface
    case = JutulCase(model, tstep, forces, state0 = state0, parameters = parameters)
    opt_config = Jutul.forces_optimization_config(model, forces, tstep, abs_min = 0.0, verbose = false)
    x0, xmin, xmax, f, g!, out = Jutul.setup_force_optimization(case, G, opt_config, verbose = false);

    der = g!(similar(x0), x0)

    @test vcat(grad_adj...) ≈ der atol = 1e-3 rtol = 1e-3
end
##
@testset "force gradients reservoir and wells" begin
    for cases in ["spe1", "egg"]
        @testset "$cases" begin
            if cases == "spe1"
                spe1_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1")
                block_backend = false
                case = setup_case_from_data_file(joinpath(spe1_dir, "SPE1.DATA"), block_backend = block_backend)[1:10]
                case = coarsen_reservoir_case(case, (3, 3, 3), setup_arg = (block_backend = block_backend,))
                injectors = [:INJ]
                producers = [:PROD]
            elseif cases == "egg"
                data_pth = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG", "EGG.DATA")
                case = setup_case_from_data_file(data_pth)
                case = case[1:10];
                case = coarsen_reservoir_case(case, (4, 4, 3), method = :ijk);
                injectors = [:INJECT1, :INJECT2, :INJECT3, :INJECT4, :INJECT5, :INJECT6, :INJECT7, :INJECT8]
                producers = [:PROD1, :PROD2, :PROD3, :PROD4]
            else
                error("Unknown case $cases")
            end
            prod_name = producers[end]
            test_force_vectorization(case.forces, case.dt, case.model)
            states, reports = simulate(case, output_substates = true, info_level = -1)
            t_tot = sum(case.dt)

            function orat_obj(model, state, dt, step_info, forces)
                orat = JutulDarcy.compute_well_qoi(model, state, forces, prod_name, SurfaceOilRateTarget)
                return dt*orat/t_tot
            end

            function prod_bhp_obj(model, state, dt, step_info, forces)
                bhp = state[prod_name].Pressure[1]
                return dt*bhp/t_tot
            end

            function cell_saturation_obj(model, state, dt, step_info, forces)
                s = state.Reservoir.Saturations[end, end]
                return dt*s/t_tot
            end

            function npv_test_obj(model, state, dt, step_info, forces)
                return JutulDarcy.npv_objective(model, state, dt, step_info, forces,
                    injectors = injectors,
                    producers = producers,
                    timesteps = case.dt
                )
            end

            function test_grad(obj; kwarg...)
                grad_eps = 1e-3
                dx = numerical_diff_forces(case.model, case.state0, case.parameters, case.forces, case.dt, obj, grad_eps; kwarg...)
                dforces, t_to_f, grad_adj = Jutul.solve_adjoint_forces(case.model, states, reports, obj, case.forces;
                    state0 = case.state0,
                    parameters = case.parameters,
                    kwarg...
                )
                for i in eachindex(dx, grad_adj)
                    @test isapprox(dx[i], grad_adj[i], atol = 1e-3, rtol = 1e-2)
                    @test norm(grad_adj, 2) ≈ norm(dx, 2) atol = 1e-3 rtol = 1e-2
                end
            end
            objectives = [cell_saturation_obj, prod_bhp_obj, npv_test_obj]
            if reservoir_model(case.model).system isa StandardBlackOilSystem
                Rs0 = sum(case.state0[:Reservoir][:Rs])

                function rs_obj(model, state, dt, step_info, forces)
                    rs = state.Reservoir.Rs
                    val = 0
                    for i in 1:length(rs)
                        val += (rs[i] - 100)^2
                    end
                    return dt*(val/(Rs0*t_tot))^2
                end

                push!(objectives, rs_obj)
            end

            for (objno, obj) in enumerate(objectives)
                @testset "Objective $obj" begin
                    test_grad(obj)
                    if objno == 1
                        # No need to do this for every combo
                        test_grad(obj, eachstep = true)
                    end
                end
            end
        end
    end
end
