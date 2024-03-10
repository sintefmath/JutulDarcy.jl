using Jutul, JutulDarcy, Test, HYPRE, MPI, PartitionedArrays

@testset "NLDD/MPI" begin
    f = JutulDarcy.GeoEnergyIO.test_input_file_path("spe1", "BENCH_SPE1.DATA")
    il = -1
    case = setup_case_from_data_file(f, split_wells = true);
    ws_base, = simulate_reservoir(case, method = :newton, info_level = il);

    ws = Dict()
    ws[:nldd], = simulate_reservoir(case, method = :nldd, info_level = il);
    ws[:nldd_4], = simulate_reservoir(case, method = :nldd, info_level = il, nldd_arg = Dict(:no_blocks => 4));

    ws[:aspen], = simulate_reservoir(case, method = :aspen, info_level = il);
    ws[:mpi], = simulate_reservoir(case, method = :newton, mode = :mpi, info_level = il);
    ws[:nldd_mpi], = simulate_reservoir(case, method = :nldd, mode = :mpi, info_level = il);

    get_prod_rate(x) = x[:PRODUCER][Symbol("Surface total rate")]
    get_grat(x) = x[:PRODUCER][Symbol("Surface gas rate")]
    get_prod_bhp(x) = x[:PRODUCER][Symbol("Bottom hole pressure")]

    for (k, ws_compare) in pairs(ws)
        @testset "$k producer" begin
            rate_ref = get_prod_rate(ws_base)
            rate = get_prod_rate(ws_compare)
            @test rate_ref ≈ rate rtol = 1e-3

            rate_ref = get_grat(ws_base)
            rate = get_grat(ws_compare)
            @test rate_ref ≈ rate rtol = 1e-3

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
        case, mrst_data = setup_case_from_mrst(name,
            split_wells = true,
            ds_max = 0.1,
            block_backend = b
        );
        push!(case.model[:Reservoir].output_variables, :Saturations)
        sim_arg = (
                    info_level = -1,
                    steps = steps,
                    block_backend = b,
                    maxstep = Inf,
                    timestepping = :iteration,
                    error_on_incomplete = true,
                    case = case, mrst_data = mrst_data
        );
        results_newton = JutulDarcy.NLDD.bench_dd(name, :fi; sim_arg...);
        its_newton = get_its(results_newton)

        function test_solver_result(result)
            @test its_newton > get_its(result)
        end
        @testset "NLDD on $name" begin
            for solve_tol_mobility in [nothing, 0.01, 0.05]
                @testset "Gauss-seidel λ_t = $solve_tol_mobility" begin
                    results_gs = JutulDarcy.NLDD.bench_dd(name, :nldd; sim_arg...,
                        gauss_seidel = true,
                        solve_tol_mobility = solve_tol_mobility);
                    test_solver_result(results_gs)
                end
                @testset "Jacobi λ_t = $solve_tol_mobility" begin
                    results_dd = JutulDarcy.NLDD.bench_dd(name, :nldd; sim_arg...,
                        gauss_seidel = false,
                        solve_tol_mobility = solve_tol_mobility);
                    test_solver_result(results_dd)
                end
            end
            if name != "spe9"
                @testset "ASPEN" begin
                    results_as = JutulDarcy.NLDD.bench_dd(name, :aspen;
                    sim_arg...);
                    test_solver_result(results_as)
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
