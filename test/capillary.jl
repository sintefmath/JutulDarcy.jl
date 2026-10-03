using Test, JutulDarcy, Jutul

@testset "Capillary pressure tables and scaling" begin
    irregular = get_1d_interpolator(
        [0.0, 0.3, 1.0], [0.0, 1.0, 2.0]; constant_dx = false)
    for constructor in (JutulDarcy.SimpleCapillaryPressure,
            JutulDarcy.ScaledCapillaryPressure)
        pc = constructor(((nothing, irregular), irregular); regions = [1, 2])
        @test typeof(pc.pc[1][1]) === typeof(irregular)
        @test pc.pc[1][1](0.4) == 0.0
        @test pc.pc[1][2](0.4) ≈ irregular(0.4)
        @test pc.pc[2][1] === irregular

        empty_pc = constructor((nothing, nothing))
        @test all(table -> table isa Jutul.LinearInterpolant,
            (empty_pc.pc[1][1], empty_pc.pc[2][1]))
    end

    with_lookup = get_1d_interpolator([0.0, 1.0], [0.0, 1.0])
    @test typeof(irregular) !== typeof(with_lookup)
    @test_logs (:warn, r"mixed types") JutulDarcy.region_wrap((irregular, with_lookup), [1, 2])

    two_phase_table = get_1d_interpolator([0.0, 1.0], [0.0, 2.0])
    saturations = reshape([0.25, 0.75], 2, 1)
    scaling = reshape([3.0], 1, 1)
    simple = JutulDarcy.SimpleCapillaryPressure((two_phase_table,))
    scaled = JutulDarcy.ScaledCapillaryPressure((two_phase_table,))
    for (reference, saturation_row) in ((1, 2), (2, 1))
        model = (system = ImmiscibleSystem((AqueousPhase(), LiquidPhase());
            reference_phase_index = reference),)
        result = zeros(1, 1)
        JutulDarcy.update_capillary_pressure!(
            result, simple, model, saturations, nothing, 1:1)
        @test result[1, 1] ≈ two_phase_table(saturations[saturation_row, 1])
        JutulDarcy.update_capillary_pressure!(
            result, scaled, model, saturations, scaling, 1:1)
        @test result[1, 1] ≈ 3*two_phase_table(saturations[saturation_row, 1])
    end

    gas_table = get_1d_interpolator([0.0, 1.0], [0.0, 4.0])
    model = (system = ImmiscibleSystem(:wog),)
    saturations = [0.2 0.6; 0.5 0.2; 0.3 0.2]
    scaling = [2.0 3.0; 4.0 5.0]
    for (constructor, scale) in (
            (JutulDarcy.SimpleCapillaryPressure, nothing),
            (JutulDarcy.ScaledCapillaryPressure, scaling))
        pc = constructor((nothing, (gas_table, gas_table)); regions = [1, 2])
        @test length(pc.pc[1]) == 2
        result = zeros(2, 2)
        JutulDarcy.update_capillary_pressure!(
            result, pc, model, saturations, scale, 1:2)
        @test all(iszero, result[1, :])
        factor = isnothing(scale) ? 1.0 : scale[2, 1]
        @test result[2, 1] ≈ factor*gas_table(saturations[3, 1])
        factor = isnothing(scale) ? 1.0 : scale[2, 2]
        @test result[2, 2] ≈ factor*gas_table(saturations[3, 2])
    end
end
