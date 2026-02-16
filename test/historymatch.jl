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
    end
end
