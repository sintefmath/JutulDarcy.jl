using Jutul, JutulDarcy, Test

function blackoil_scaling_test_model()
    pressure = [1.0e7, 2.0e7, 3.0e7]
    table(values) = get_1d_interpolator(pressure, values; cap_end = false)
    system = StandardBlackOilSystem(
        rs_max = (table([10.0, 100.0, 80.0]), table([20.0, 200.0, 160.0])),
        rv_max = (table([0.001, 0.01, 0.008]), table([0.002, 0.02, 0.016])))
    domain = reservoir_domain(CartesianMesh((3,), (3.0,)))
    return SimulationModel(domain, system)
end

function blackoil_scaling_test_state()
    return (
        Pressure = fill(2.0e7, 3),
        ImmiscibleSaturation = fill(0.1, 3),
        BlackOilUnknown = [
            BlackOilX(20.0, JutulDarcy.OilOnly),
            BlackOilX(0.002, JutulDarcy.GasOnly),
            BlackOilX(0.3, JutulDarcy.OilAndGas)
        ]
    )
end

@testset "Black-oil variable scaling" begin
    model = blackoil_scaling_test_model()
    variable = model.primary_variables[:BlackOilUnknown]
    @test variable.rs_scale == 200.0
    @test variable.rv_scale == 0.02
    @test BlackOilUnknown().rs_scale == BlackOilUnknown().rv_scale == 1.0
    @test BlackOilUnknown(StandardBlackOilSystem()).rs_scale == 1.0
    @test BlackOilUnknown(StandardBlackOilSystem()).rv_scale == 1.0
    @test BlackOilUnknown(StandardBlackOilSystem(rs_max = model.system.rs_max)).rv_scale == 1.0
    @test BlackOilUnknown(StandardBlackOilSystem(rv_max = model.system.rv_max)).rs_scale == 1.0
    custom = BlackOilUnknown(model.system; rs_scale = 50.0, rv_scale = 0.005)
    @test custom.rs_scale == 50.0
    @test custom.rv_scale == 0.005
    @test BlackOilUnknown(rs_scale = 200.0f0) isa BlackOilUnknown{Float64}
    @test BlackOilUnknown{Float32}().rs_scale === 1.0f0
    @test BlackOilUnknown(Inf, 0.2).rv_scale == 1.0
    callable_system = StandardBlackOilSystem(rs_max = p -> 100.0;
        eps_rs = 1e-8)
    @test BlackOilUnknown(callable_system).rs_scale == 1.0
    setup_model = setup_reservoir_model(model.data_domain, model.system;
        ds_max = 0.15, dr_max = 0.3)
    setup_variable = reservoir_model(setup_model).primary_variables[:BlackOilUnknown]
    @test setup_variable.rs_scale == variable.rs_scale
    @test setup_variable.rv_scale == variable.rv_scale
    @test setup_variable.ds_max == 0.15
    @test setup_variable.dr_max == 0.3

    @testset "Physical increments and limits" begin
        for w in (1.0, 0.5)
            state = blackoil_scaling_test_state()
            Jutul.update_primary_variable!(state, variable, :BlackOilUnknown,
                model, fill(0.1, 3), w)
            @test [x.val for x in state.BlackOilUnknown] ≈
                [20.0 + 20.0*w, 0.002 + 0.002*w, 0.3 + 0.1*w]
        end
        limited = BlackOilUnknown(model.system; dr_max = 0.1, ds_max = 0.05)
        state = blackoil_scaling_test_state()
        Jutul.update_primary_variable!(state, limited, :BlackOilUnknown,
            model, ones(3), 1.0)
        @test [x.val for x in state.BlackOilUnknown] ≈ [30.0, 0.003, 0.35]

        state = blackoil_scaling_test_state()
        Jutul.update_primary_variable!(state, variable, :BlackOilUnknown,
            model, fill(-1.0, 3), 1.0)
        @test state.BlackOilUnknown[1].val == 0.0
        @test state.BlackOilUnknown[2].val == 0.0
    end

    @testset "Derivatives follow phase switching" begin
        state = Dict(pairs(blackoil_scaling_test_state()))
        Jutul.initialize_primary_variable_ad!(state, model, variable,
            :BlackOilUnknown, 3; offset = 1, tag = Jutul.Cells())
        state = (; pairs(state)...)
        unknown = state[:BlackOilUnknown]
        ratios = similar(unknown, typeof(unknown[1].val))
        function check_ratios()
            JutulDarcy.update_rs!(ratios, model.secondary_variables[:Rs],
                model, state[:Pressure], unknown, 1:3)
            for i in 1:3
                expected = unknown[i].phases_present == JutulDarcy.OilOnly ? variable.rs_scale : 0.0
                @test Jutul.ForwardDiff.partials(ratios[i])[2] == expected
                if expected > 0
                    @test Jutul.value(ratios[i]) == Jutul.value(unknown[i].val)
                end
            end
            JutulDarcy.update_rv!(ratios, model.secondary_variables[:Rv],
                model, state[:Pressure], unknown, 1:3)
            for i in 1:3
                expected = unknown[i].phases_present == JutulDarcy.GasOnly ? variable.rv_scale : 0.0
                @test Jutul.ForwardDiff.partials(ratios[i])[2] == expected
                if expected > 0
                    @test Jutul.value(ratios[i]) == Jutul.value(unknown[i].val)
                end
            end
        end
        check_ratios()
        # Both undersaturated variables cross their saturation limits.
        Jutul.update_primary_variable!(state, variable, :BlackOilUnknown,
            model, [1.0, 1.0, -1.0], 1.0)
        @test unknown[1].phases_present == JutulDarcy.OilAndGas
        @test unknown[2].phases_present == JutulDarcy.OilAndGas
        check_ratios()
        # The two-phase cells now lose gas and oil, respectively.
        Jutul.update_primary_variable!(state, variable, :BlackOilUnknown,
            model, [-1.0, 1.0, 0.0], 1.0)
        @test unknown[1].phases_present == JutulDarcy.OilOnly
        @test unknown[2].phases_present == JutulDarcy.GasOnly
        check_ratios()
    end

    @testset "ASPEN increments" begin
        state = blackoil_scaling_test_state()
        buffer = zeros(3)
        extra = zeros(3, 3)
        JutulDarcy.NLDD.transfer_pvar!(buffer, model, state, 1:3, false,
            state.BlackOilUnknown, variable, extra, 1.0, 0)
        Jutul.update_primary_variable!(state, variable, :BlackOilUnknown,
            model, fill(0.1, 3), 1.0)
        JutulDarcy.NLDD.transfer_pvar!(buffer, model, state, 1:3, true,
            state.BlackOilUnknown, variable, extra, 1.0, 0)
        @test buffer ≈ fill(-0.1, 3)
    end

    @testset "Increment reports and backend adaptation" begin
        state = blackoil_scaling_test_state()
        for increments in (fill(0.1, 3), fill(0.1, 3)')
            report = Jutul.increment_norm(increments, state, model,
                state.BlackOilUnknown, variable)
            @test report.rs_sum == report.rs_max ≈ 20.0
            @test report.rv_sum == report.rv_max ≈ 0.002
            @test report.sg_sum == report.sg_max ≈ 0.1
        end
        context = KernelAbstractionsContext(
            Jutul.KernelExecution.KernelAbstractions.CPU(); float_type = Float32)
        adapted = JutulDarcy.Adapt.adapt(context, variable)
        @test adapted isa BlackOilUnknown{Float32}
        @test adapted.rs_scale === 200.0f0
        @test adapted.rv_scale === 0.02f0
    end
end
