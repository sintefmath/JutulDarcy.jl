using Jutul, JutulDarcy
using JutulDarcy.HistoryMatching

using Test
@testset "HistoryMatching Interface" begin
    spe1_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1")
    case_spe1 = setup_case_from_data_file(joinpath(spe1_dir, "SPE1.DATA"))
    res_spe1 = simulate_reservoir(case_spe1, info_level = -1)

    yr = si_unit(:year)
    hm = HistoryMatching.HistoryMatch(case_spe1)
    @test ismissing(hm.periods)
    hm = HistoryMatching.HistoryMatch(case_spe1, periods = [5.0, 7.5].*yr)
    @test length(hm.periods) == 2
    @test hm.periods[1].start == 0.0
    @test hm.periods[1].stop == 5.0*yr
    @test hm.periods[2].start == 5.0*yr
    @test hm.periods[2].stop == 7.5*yr

    hm = HistoryMatching.HistoryMatch(case_spe1, periods = [(0, 5.0yr), (7.5yr, 8.0yr)], period_weights = [0.3, 0.7])
    for i in findall(x -> x < 5yr, cumsum(case_spe1.dt))
        @test i in hm.periods[1].step_idx
    end

    for i in findall(x -> x > 7.5yr && x < 8yr, cumsum(case_spe1.dt))
        @test i in hm.periods[2].step_idx
    end
    @test length(hm.periods) == 2
    @test hm.periods[1].start == 0.0
    @test hm.periods[1].stop == 5.0*yr
    @test hm.periods[2].start == 7.5*yr
    @test hm.periods[2].stop == 8.0*yr
    @test hm.periods[1].weights[hm.periods[1].step_idx[1]] <= 0.3
    @test hm.periods[1].weights[hm.periods[1].step_idx[2]] == 0.3
    @test hm.periods[2].weights[hm.periods[2].step_idx[1]] <= 0.7
    @test hm.periods[2].weights[hm.periods[2].step_idx[2]] == 0.7

    hm = HistoryMatching.HistoryMatch(case_spe1, res_spe1)

    @test_throws "" HistoryMatching.match_injectors!(hm, :badval)

    HistoryMatching.match_injectors!(hm, :bhp)
    HistoryMatching.match_producers!(hm, :grat, weight = 2.0)

    @test HistoryMatching.get_injectors(hm) == [:INJ]
    @test HistoryMatching.get_producers(hm) == [:PROD]

    @testset "GlobalObjective" begin
        obj = history_match_objective(case_spe1, res_spe1, is_global = true)
        match_injectors!(obj, :bhp)
        match_producers!(obj, :grat, weight = 2.0)
        @test evaluate_match(obj, res_spe1) ≈ 0.0 atol = 1e-12
        ##
        match_producers!(obj, :cumulative_water, weight = 2.0)
        @test evaluate_match(obj, res_spe1) ≈ 0.0 atol = 1e-12
    end
    ##
    @testset "SumObjective" begin
        @test_throws "" obj_local = history_match_objective(case_spe1, res_spe1, is_global = false, periods = [5.0, 7.5].*yr)
        obj_local = history_match_objective(case_spe1, res_spe1, is_global = false)
        match_injectors!(obj_local, :bhp)
        match_producers!(obj_local, :grat, weight = 2.0)
        @test_throws "Cumulative well matches are not supported for SumHistoryMatchObjective." match_producers!(obj_local, :cumulative_water, weight = 2.0)
        @test evaluate_match(obj_local, res_spe1) ≈ 0.0 atol = 1e-12
        nstep = length(case_spe1.dt)
        obj2 = history_match_objective(case_spe1, res_spe1, is_global = false)
        match_injectors!(obj2, :bhp, weight = rand(nstep))
        # Mismatching weight vector length should throw an error
        @test_throws "Weight vector length must match number of simulation report steps in case" match_injectors!(obj2, :bhp, weight = rand(nstep-1))
    end

    @testset "WWCT/WWGR/WGOR/WGLR" begin
        spe1_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1")
        spe1_file = joinpath(spe1_dir, "SPE1.DATA")
        spe = JutulDarcy.GeoEnergyIO.parse_data_file(spe1_file)
        @. spe["PROPS"]["SWOF"][1][:, 2] *= 1e5
        case = setup_case_from_data_file(spe)

        res = simulate_reservoir(case, info_level = -1)
        res2 = deepcopy(res)
        res1 = deepcopy(res)
        s = res1.summary
        s["VALUES"]["WELLS"]["PROD"]["WWCT"] .= 1.0
        s["VALUES"]["WELLS"]["PROD"]["WWGR"] .= 1.0
        s["VALUES"]["WELLS"]["PROD"]["WGOR"] .= 1.0
        s["VALUES"]["WELLS"]["PROD"]["WGLR"] .= 1.0

        for glob in [true, false]
            for wkey in ["WWCT", "WWGR", "WGOR", "WGLR"]
                obj1 = history_match_objective(case, res1, is_global = glob, scale = 1.0)
                match_producers!(obj1, wkey, scale = 1.0)

                obj2 = history_match_objective(case, res2, is_global = glob, scale = 1.0)
                match_producers!(obj2, wkey, scale = 1.0)

                obs = res1.summary["VALUES"]["WELLS"]["PROD"][wkey]
                ref = res2.summary["VALUES"]["WELLS"]["PROD"][wkey]
                w = case.dt
                delta = (abs.(obs .- ref).^2).*w
                perturbed_match = evaluate_match(obj1, res1)
                ref_value = sum(delta)
                @test perturbed_match ≈ ref_value rtol = 1e-3
                zero_match = evaluate_match(obj2, res2)
                @test zero_match ≈ 0.0 atol = 1e-12
            end
        end
    end
end
##
