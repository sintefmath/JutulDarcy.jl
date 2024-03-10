using Test

function compare_states(states, states_other, rtol)
    @test length(states) == length(states_other)
    for i in eachindex(states)
        @testset "step $i" begin
            if length(states_other) < i
                @test false
            else
                for (k, v) in states_other[i]
                    if k == :Pressure
                        # Base tests might be singular and pressure may not
                        # be unique. We anyway test TotalMasses so this should
                        # be fine to skip.
                        continue
                    end
                    ref = states[i][k]
                    @test norm(ref - v)/norm(ref) < rtol
                end
            end
        end
    end
end

function test_aspen_simple(method = :aspen, mesh = "cart_2x1", block = false, casename = "two_phase_simple"; info_level = -1, dt = repeat([10.0], 10), max_timestep_cuts = 0, extra_arg...)
    if block
        context = DefaultContext(matrix_layout = BlockMajorLayout())
        lsolve = GenericKrylov(:bicgstab, relative_tolerance = 1e-12, preconditioner = CPRPreconditioner())
    else
        context = DefaultContext()
        lsolve = nothing
    end
    if mesh == "cart_100x1"
        grid = CartesianMesh((100, 1), (10000.0, 1.0))
        nc = number_of_cells(grid)
        N = 5
        p = partition_uniform_1d(nc, N)
    elseif mesh == "cart_3x4"
        grid = CartesianMesh((3, 4), (1000.0, 1000.0))
        nc = number_of_cells(grid)
        N = 5
        p = partition_uniform_1d(nc, N)
    elseif mesh == "cart_2x2"
        grid = CartesianMesh((2, 2), (1000.0, 1000.0))
        nc = number_of_cells(grid)
        N = 2
        p = partition_uniform_1d(nc, N)
    elseif mesh == "cart_2x2_trivial_partition"
        grid = CartesianMesh((2, 2), (1000.0, 1000.0))
        nc = number_of_cells(grid)
        N = 4
        p = partition_uniform_1d(nc, N)
    elseif mesh == "cart_1x2_trivial_partition"
        grid = CartesianMesh((1, 2), (1000.0, 1000.0))
        nc = number_of_cells(grid)
        N = 2
        p = [2, 1]
    elseif mesh == "cart_2x2_oneblock"
        grid = CartesianMesh((2, 2), (1000.0, 1000.0))
        nc = number_of_cells(grid)
        N = 1
        p = ones(Int64, nc)
    elseif mesh == "cart_6x1"
        grid = CartesianMesh((6, 1), (10000.0, 1.0))
        p = [1, 1, 3, 3, 2, 2]
    elseif mesh == "cart_4x1"
        grid = CartesianMesh((4, 1), (10000.0, 1.0))
        p = [1, 1, 2, 2]
    elseif mesh == "cart_2x1"
        grid = CartesianMesh((2, 1), (2.0, 1.0))
        p = [2, 1]
    else
        error("Unsupported mesh $mesh")
    end
    part = SimplePartition(p)
    state0, model, prm, f, t = get_test_setup(grid, case_name = casename, pvfrac = 100.0, timesteps = dt, context = context)
    sim = Simulator(model, state0 = state0, parameters = prm)

    make_config(sim; kwarg...) = simulator_config(sim; kwarg..., linear_solver = lsolve,
                                                             max_timestep_cuts = max_timestep_cuts,
                                                             info_level = info_level, extra_arg...)
    # Simulate Newton
    cfg = make_config(sim)
    states, reports = simulate(sim, t, forces = f, config = cfg);
    # Simulate NLDD/ASPEN
    sim_dd = NLDDSimulator(model, part, state0 = state0, parameters = prm);
    cfg_dd = make_config(sim_dd, method = method)
    submodels = map(x -> x.model, sim_dd.subdomain_simulators)
    # f_dd = aspen_forces(f, submodels)
    states_dd, reports_dd = simulate(sim_dd, t, forces = f, config = cfg_dd);

    rtol = 0.01
    @testset "$method" begin
        compare_states(states, states_dd, rtol)
    end
    return (states_dd, states, reports_dd, reports, sim, sim_dd)
end

function test_aspen_physics(physics::Vector)
    meshes = ["cart_100x1", "cart_3x4", "cart_6x1", "cart_4x1", "cart_2x1"]
    @testset "ASPEN/NLDD all physics" begin
        for p in physics
            for m in meshes
                for method in [:aspen, :nldd]
                    @testset "$method - $m $p" begin
                        test_aspen_simple(method, m, false, p)
                        # test_aspen_simple(method, m, true, p)
                    end
                end
            end
        end
    end
end

test_aspen_physics(s::String) = test_aspen_physics([s])

function test_aspen_physics()
    physics = [
        # "single_phase_simple",
        "two_phase_simple",
        "simple_compositional_fake_wells",
        "two_phase_fake_wells",
        "three_phase_fake_wells"
    ]
    test_aspen_physics(physics)
end
