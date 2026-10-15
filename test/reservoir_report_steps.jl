using Test, JutulDarcy, Jutul

@testset "reservoir_report_steps" begin
    case = JutulDarcy.setup_mini_wellcase(Val(:bo_spe1), nstep = 10)
    nstep = length(case.dt)
    steps = [2, 5, 7]
    ref = simulate_reservoir(case, info_level = -1)

    for pth in (nothing, tempname())
        res = simulate_reservoir(case,
            info_level = -1,
            output_path = pth,
            reservoir_report_steps = steps
        )
        @testset "output_path = $pth" begin
            @test length(res.states) == nstep
            for i in 1:nstep
                if i in steps
                    @test res.states[i] == ref.states[i]
                else
                    @test ismissing(res.states[i])
                end
                # Well states are always stored
                @test haskey(res.result.states[i], :Producer)
                @test haskey(res.result.states[i], :Facility)
            end
            for (k, v) in pairs(ref.wells.wells)
                for (wk, wv) in pairs(v)
                    @test res.wells[k][wk] == wv
                end
            end
            field = res.summary["VALUES"]["FIELD"]
            field_ref = ref.summary["VALUES"]["FIELD"]
            @test field["FOPR"] == field_ref["FOPR"]
            fpr = field["FPR"]
            @test fpr[steps] == field_ref["FPR"][steps]
            @test all(isnan, fpr[setdiff(1:nstep, steps)])
            @test occursin("stored for 3 states", sprint(show, MIME("text/plain"), res))
        end
    end

    @testset "Composes with output_function" begin
        function F(state, report)
            state = copy(state)
            state[:Extra] = true
            return state
        end
        res = simulate_reservoir(case,
            info_level = -1,
            output_function = F,
            reservoir_report_steps = steps
        )
        @test all(s -> haskey(s, :Extra), res.result.states)
        @test count(!ismissing, res.states) == length(steps)
    end

    @testset "Restart" begin
        pth = tempname()
        # Simulate first part of schedule
        case_part = case[1:6]
        simulate_reservoir(case_part,
            info_level = -1,
            output_path = pth,
            reservoir_report_steps = steps
        )
        err = try
            simulate_reservoir(case, info_level = -1, output_path = pth, restart = true)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("Use `restart = 6`", err.msg)
        @test_throws ErrorException simulate_reservoir(case, info_level = -1, output_path = pth, restart = 4)
        res = simulate_reservoir(case, info_level = -1, output_path = pth, restart = 6)
        @test length(res.states) == nstep
        @test res.states[7][:Pressure] ≈ ref.states[7][:Pressure]
    end
end
