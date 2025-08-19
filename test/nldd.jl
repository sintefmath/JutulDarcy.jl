using Jutul, JutulDarcy, Test, HYPRE, MPI, PartitionedArrays

@testset "NLDD/MPI" begin
    f = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
    il = -1
    case = setup_case_from_data_file(f, split_wells = true);
    ws_base, = simulate_reservoir(case, method = :newton, info_level = il);

    ws = Dict()
    ws[:nldd], = simulate_reservoir(case, method = :nldd, info_level = il);
    ws[:nldd_4], = simulate_reservoir(case, method = :nldd, info_level = il, nldd_arg = Dict(:no_blocks => 4));

    ws[:aspen], = simulate_reservoir(case, method = :aspen, info_level = il);
    ws[:mpi], = simulate_reservoir(case, method = :newton, mode = :mpi, info_level = il);
    ws[:nldd_mpi], = simulate_reservoir(case, method = :nldd, mode = :mpi, info_level = il);

    get_prod_rate(x) = x[:PROD, :rate]
    get_grat(x) = x[:PROD, :grat]
    get_prod_bhp(x) = x[:PROD, :bhp]

    for (k, ws_compare) in pairs(ws)
        @testset "$k producer" begin
            rate_ref = get_prod_rate(ws_base)
            rate = get_prod_rate(ws_compare)
            @test rate_ref ≈ rate rtol = 1e-2

            rate_ref = get_grat(ws_base)
            rate = get_grat(ws_compare)
            @test rate_ref ≈ rate rtol = 1e-2

            pbhp_ref = get_prod_bhp(ws_base)
            pbhp = get_prod_bhp(ws_compare)
            @test pbhp_ref ≈ pbhp rtol = 1e-2
        end
    end
end
@testset "NLDD/ASPEN integration tests" begin
    JutulDarcy.NLDD.test_aspen_physics();
    ##
    function test_on_mrst_case(name, b = true; steps = :full)

        function get_its(result)
            Jutul.report_stats(result.result.reports).newtons
        end
        ds_max = 0.1
        reports_dd = states = states_as = states_gs = []
        case, = setup_case_from_mrst(name,
            split_wells = true,
            ds_max = 0.1,
            block_backend = b
        );
        if steps != :full
            case = case[steps]
        end
        push!(case.model[:Reservoir].output_variables, :Saturations)
        sim_arg = (
            info_level = -1,
            max_timestep = Inf,
            timesteps = :iteration,
            error_on_incomplete = true
        );
        results_newton = simulate_reservoir(case, method = :newton; sim_arg...);
        its_newton = get_its(results_newton)

        function test_solver_result(result)
            @test its_newton >= get_its(result)
        end
        @testset "NLDD on $name" begin
            for solve_tol_mobility in [nothing, 0.01, 0.05]
                @testset "Gauss-seidel λ_t = $solve_tol_mobility" begin
                    results_gs = simulate_reservoir(case, method = :nldd; sim_arg...,
                        gauss_seidel = true,
                        solve_tol_mobility = solve_tol_mobility);
                    test_solver_result(results_gs)
                end
                @testset "Jacobi λ_t = $solve_tol_mobility" begin
                    results_dd = simulate_reservoir(case, method = :nldd; sim_arg...,
                        gauss_seidel = false,
                        solve_tol_mobility = solve_tol_mobility);
                    test_solver_result(results_dd)
                end
            end
            if name != "spe9"
                @testset "ASPEN" begin
                    results_as = simulate_reservoir(case, method = :aspen;
                    sim_arg...);
                    test_solver_result(results_as)
                end
            end
        end
    end

    @testset "NLDD with BCs" begin
        nc = 4
        function setup_bc_case(; thermal = false)
            # Model
            grid = CartesianMesh((nc, 1, nc), (10.0, 1.0, 10.0))
            domain = reservoir_domain(grid)
            if thermal
                sys = :geothermal
            else
                sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
            end
            model, = setup_reservoir_model(domain, sys)
            # Initial state and BCs
            p0, T0 = 10.0si_unit(:bar), convert_to_si(20.0, :Celsius)
            state0 = setup_reservoir_state(model; Pressure = p0, Temperature = T0)
            bc = flow_boundary_condition([1, nc^2], domain,
                [2*p0, 0.5*p0], [T0 + 75.0, T0])
            forces = setup_reservoir_forces(model; bc = bc)
            # Case
            case = JutulCase(model, [1si_unit(:day)], forces, state0 = state0)
            return case
        end

        function get_its(result)
            Jutul.report_stats(result.result.reports).newtons
        end

        args = (info_level = -1, max_timestep = Inf)
        # for thermal = [false, true] # Fails for thermal = false (single-phase w/BCs)
        for thermal = [true]
            case = setup_bc_case(thermal = thermal)
            results_newton = simulate_reservoir(case;
            method = :newton, args...)
            for nb = [1,4]
                results_nldd = simulate_reservoir(case;
                    method = :nldd, nldd_arg = (no_blocks = nb, ), args...)
                if nb == 1
                    # Simulation should converge in a single global iteration
                    # with a single NLDD partition
                    @test get_its(results_nldd) == 1
                end
                for k = [:Pressure, :Temperature]
                    # Compare NLDD to Newton
                    v_newton = results_newton.states[end][k]
                    v_nldd = results_nldd.states[end][k]
                    @test all(isapprox.(v_nldd, v_newton, rtol=1e-3))
                end
            end
        end

    end

    @testset "NLDD on MRST cases" begin
        for use_blocks in [true, false]
            test_on_mrst_case("spe1", use_blocks)
        end
        test_on_mrst_case("egg", true, steps = 1:5)
        test_on_mrst_case("spe9", true, steps = 1:5)
    end
end
