using JutulDarcy, Jutul, Test, LinearAlgebra
import Jutul.DictOptimization: finite_difference_gradient_entry

@testset "Inner adjoint tests" begin
    function test_physics(phys;
            block_backend = false,
            case_broken = false,
            di_broken = false,
            atol = 1e-8,
            rtol = 1e-2,
            lumping = missing,
            add_let = false
        )
        tmp = JutulDarcy.simulate_mini_wellcase(phys, block_backend = block_backend);
        s = tmp.setup;
        m = s[:model]
        prm = s[:parameters]
        if add_let
            rmodel = reservoir_model(m)
            kr = JutulDarcy.ParametricLETRelativePermeabilities()
            set_secondary_variables!(rmodel, RelativePermeabilities = kr)
            JutulDarcy.add_relperm_parameters!(rmodel)
            prm = setup_parameters(m)
        end
        case = JutulCase(m, s[:dt], s[:forces], state0 = s[:state0], parameters = prm)
        _, states = simulate_reservoir(case, info_level = -1)

        step_times = cumsum(case.dt)
        total_time = step_times[end]
        function pdiff(p, p0)
            v = 0.0
            for i in eachindex(p)
                v += (p[i] - p0[i])^2
            end
            return v
        end

        function mismatch_objective_p(m, s, dt, step_info, forces)
            t = step_info[:time] + dt
            step = findmin(x -> abs(x - t), step_times)[2]
            p = s[:Reservoir][:Pressure]
            v = pdiff(p, states[step][:Pressure])
            return (dt/total_time)*(v/(si_unit(:bar)*100)^2)
        end

        function setup_function(d, step_info = missing)
            new_case = deepcopy(case)
            T = new_case.parameters[:Reservoir][:Transmissibilities]
            p0 = new_case.state0[:Reservoir][:Pressure]
            new_case.state0[:Reservoir][:Pressure] = p0.*d["pmult"]
            new_case.parameters[:Reservoir][:Transmissibilities] = T.*d["transmult"]
            if add_let
                WettingLET = new_case.parameters[:Reservoir][:WettingLET]
                NonWettingLET = new_case.parameters[:Reservoir][:NonWettingLET]
                new_case.parameters[:Reservoir][:WettingLET] = WettingLET .* d["WettingLET"]
                new_case.parameters[:Reservoir][:NonWettingLET] = NonWettingLET .* d["NonWettingLET"]
            end
            return new_case
        end

        rmodel = reservoir_model(case)
        nf = number_of_faces(rmodel.domain)
        nc = number_of_cells(rmodel.domain)
        prm = Dict{String, Any}(
            "transmult" => 1.0 .+ 0.5*rand(nf),
            "pmult" => fill(1.5, nc)
        )
        if add_let
            prm["WettingLET"] = fill(1.0, 3, nc)
            prm["NonWettingLET"] = fill(1.0, 3, nc)
        end

        dprm = setup_reservoir_dict_optimization(prm, setup_function, verbose = false)
        free_optimization_parameter!(dprm, "transmult", abs_max = 5.0, abs_min = 0.1)
        free_optimization_parameter!(dprm, "pmult", abs_max = 5.0, abs_min = 0.1, lumping = lumping)

        if add_let
            if ismissing(lumping)
                lumping_let = missing
            else
                lumping_let = repeat(lumping', 3, 1)
            end
            free_optimization_parameter!(dprm, "WettingLET", abs_max = 5.0, abs_min = 0.1, lumping = lumping_let)
            free_optimization_parameter!(dprm, "NonWettingLET", abs_max = 5.0, abs_min = 0.1)
        end

        P = Jutul.DictOptimization.JutulOptimizationProblem(dprm, mismatch_objective_p, setup_function, info_level = -1, deps = :parameters_and_state0, deps_ad = :jutul)
        P_case = Jutul.DictOptimization.JutulOptimizationProblem(dprm, mismatch_objective_p, setup_function, info_level = -1, deps = :case)
        P_di = Jutul.DictOptimization.JutulOptimizationProblem(dprm, mismatch_objective_p, setup_function, info_level = -1, deps = :parameters_and_state0, deps_ad = :di)
        x0 = copy(P.x0)

        f, dfdx_jutul = P(x0)
        f, dfdx_di = P_di(x0)
        f, dfdx_case = P_case(x0)
        ##
        dfdx_fd = similar(dfdx_jutul)
        for index in eachindex(dfdx_jutul)
            dfdx_fd[index] = finite_difference_gradient_entry(P, x0, index = index)
        end
        function test_ok(val, broken)
            @test length(val) == length(dfdx_fd)
            normdiff = norm(val-dfdx_fd, 2)/norm(dfdx_fd, 2)
            if broken
                @test_broken normdiff < 1e-3
                @test_broken all(isapprox.(dfdx_fd, val, rtol = rtol, atol = atol))
            else
                @test normdiff < 1e-3
                @test all(isapprox.(dfdx_fd, val, rtol = rtol, atol = atol))
            end
        end
        @testset "Jutul inner adjoint" begin
            test_ok(dfdx_jutul, false)
        end
        @testset "DI inner adjoint" begin
            test_ok(dfdx_di, di_broken)
        end
        @testset "DI outer adjoint (case)" begin
            test_ok(dfdx_case, case_broken)
        end
    end
    phys = :immiscible_2ph
    @testset "$phys" begin
        test_physics(phys; block_backend = true)
        test_physics(phys; block_backend = false)

        @testset "Lumping" begin
            test_physics(phys; block_backend = true, lumping = [2, 1, 1])
        end
        @testset "LET" begin
            test_physics(phys; block_backend = true, lumping = [2, 1, 1], add_let = true, di_broken = false)
        end
    end
    phys = :compositional_2ph_3c
    @testset "$phys" begin
        test_physics(phys; block_backend = true)
        test_physics(phys; block_backend = false)
    end
    phys = :bo_spe1
    @testset "$phys" begin
        test_physics(phys; block_backend = true, case_broken = false, di_broken = true)
        test_physics(phys; block_backend = false, case_broken = false, di_broken = true)
    end
    phys = :single_phase
    @testset "$phys" begin
        test_physics(phys; block_backend = false, case_broken = false, di_broken = false)
    end
end;
