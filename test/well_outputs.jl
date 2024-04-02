using Test, JutulDarcy

@testset "Controls to symbols" begin
    @test JutulDarcy.translate_target_to_symbol(BottomHolePressureTarget(1.0)) == :bhp
    @test JutulDarcy.translate_target_to_symbol(TotalRateTarget(1.0)) == :rate
    @test JutulDarcy.translate_target_to_symbol(SurfaceWaterRateTarget(1.0)) == :wrat
    @test JutulDarcy.translate_target_to_symbol(SurfaceLiquidRateTarget(1.0)) == :lrat
    @test JutulDarcy.translate_target_to_symbol(SurfaceOilRateTarget(1.0)) == :orat
    @test JutulDarcy.translate_target_to_symbol(SurfaceGasRateTarget(1.0)) == :grat
    @test JutulDarcy.translate_target_to_symbol(ReservoirVoidageTarget(1.0, (0.2, 0.5, 0.3))) == :resv
    @test JutulDarcy.translate_target_to_symbol(HistoricalReservoirVoidageTarget(1.0, (0.2, 0.5, 0.3))) == :resv_history
end
