using Test, Jutul, JutulDarcy
##
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
            @testset "ASPEN" begin
                results_as = JutulDarcy.NLDD.bench_dd(name, :aspen;
                sim_arg...);
                test_solver_result(results_as)
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
