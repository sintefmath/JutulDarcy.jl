using Test, JutulDarcy, Jutul, Dates

@testset "Controls to symbols" begin
    @test JutulDarcy.translate_target_to_symbol(BottomHolePressureTarget(1.0)) == :bhp
    @test JutulDarcy.translate_target_to_symbol(TotalRateTarget(1.0)) == :rate
    @test JutulDarcy.translate_target_to_symbol(SurfaceWaterRateTarget(1.0)) == :wrat
    @test JutulDarcy.translate_target_to_symbol(SurfaceLiquidRateTarget(1.0)) == :lrat
    @test JutulDarcy.translate_target_to_symbol(SurfaceOilRateTarget(1.0)) == :orat
    @test JutulDarcy.translate_target_to_symbol(SurfaceGasRateTarget(1.0)) == :grat
    @test JutulDarcy.translate_target_to_symbol(ReservoirVoidageTarget(1.0, (0.2, 0.5, 0.3))) == :resv
    @test JutulDarcy.translate_target_to_symbol(HistoricalReservoirVoidageTarget(1.0, (0.2, 0.5, 0.3))) == :resv_history
    @test JutulDarcy.translate_target_to_symbol(TotalMassRateTarget(1.0)) == :mrat
end

@testset "WellResults" begin
    n = 10
    t = collect(range(1.0, 10.0, n))
    W1 = Dict{Symbol, Vector}(:rate => rand(n), :bhp => 1e8.*rand(n), :unknown => -rand(n))
    W2 = Dict{Symbol, Vector}(:rate => rand(n), :bhp => 1e8.*rand(n), :unknown => -rand(n))

    ws = Dict(
        :W1 => W1,
        :W2 => W2
    )

    res = JutulDarcy.WellResults(t, ws)

    @test res[:W1] == W1
    @test res[:W2] == W2
    @test res[:W2, :bhp] == W2[:bhp]
    @test res[:W2, [:bhp, :rate]] == res[:W2, (:bhp, :rate)] == hcat(W2[:bhp], W2[:rate])
    @test res[[:W1, :W2], :bhp] == res[(:W1, :W2), :bhp] == hcat(W1[:bhp], W2[:bhp])
    @test res[[:W1, :W2], [:bhp, :rate]] == [hcat(W1[:bhp], W1[:rate]), hcat(W2[:bhp], W2[:rate])]

    res_date = JutulDarcy.WellResults(t, ws, today())
    JutulDarcy.well_result_dates(res)
    JutulDarcy.well_result_dates(res_date)

    @test JutulDarcy.print_well_result_table(res, :W1) == nothing
    @test JutulDarcy.print_well_result_table(res_date, :W1) == nothing

    @test JutulDarcy.print_well_result_table(res_date, (:W1, :W2), :bhp) == nothing
    @test JutulDarcy.print_well_result_table(res_date, (:W1, :W2)) == nothing
    @test res_date(:W1) == nothing
end
