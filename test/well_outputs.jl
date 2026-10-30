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

@testset "Rate targets respect their lower pressure constraint" begin
    day, bar = si_units(:day, :bar)
    ctrl = ProducerControl(SurfaceWaterRateTarget(-125/day))
    cfg = JutulDarcy.WellGroupConfiguration([:P], Dict(:P => ctrl), Dict(:P => (bhp = bar,)))
    state = (WellGroupConfiguration = cfg,)
    function control_residual(q, p)
        cond = (surface_aqueous_rate = q, bottom_hole_pressure = p)
        return JutulDarcy.well_control_equation(ctrl, cond, :P, nothing, state)
    end
    @test control_residual(-125/day, 180bar) ≈ 0 atol = 1e-12
    @test control_residual(0.0, bar) == 0
    @test control_residual(0.0, 180bar) > 0
    @test control_residual(-126/day, 180bar) < 0
    @test control_residual(0.0, 0.5bar) < 0
    # The pressure constraint supplies a Jacobian entry even with no water.
    @test Jutul.ForwardDiff.derivative(p -> control_residual(0.0, p), 180bar) != 0
    @test all(isfinite, Jutul.ForwardDiff.gradient(v -> control_residual(v...), [-125/day, bar]))
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

    @test isnothing(JutulDarcy.print_well_result_table(res, :W1))
    @test isnothing(JutulDarcy.print_well_result_table(res_date, :W1))

    @test isnothing(JutulDarcy.print_well_result_table(res_date, (:W1, :W2), :bhp))
    @test isnothing(JutulDarcy.print_well_result_table(res_date, (:W1, :W2)))
    @test isnothing(res_date(:W1))
end
