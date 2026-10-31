using Jutul, JutulDarcy, Test

function limit_test_cond(bhp, qmass)
    return JutulDarcy.FacilityVariablesForWell(1, :Well, bhp, qmass, 0.0, 0.0, 0.0)
end

@testset "Mass rate limits" begin
    @testset "Producer" begin
        limits = (bhp = 50e5, mrat = 100.0)
        # Mass rate target, BHP limit violated
        ctrl = ProducerControl(TotalMassRateTarget(-100.0))
        new_ctrl, changed = JutulDarcy.check_well_limits(limits, limit_test_cond(40e5, -100.0), ctrl)
        @test changed
        @test new_ctrl.target isa BottomHolePressureTarget
        @test new_ctrl.target.value == 50e5
        # Mass rate target, BHP limit satisfied
        new_ctrl, changed = JutulDarcy.check_well_limits(limits, limit_test_cond(60e5, -100.0), ctrl)
        @test !changed
        # BHP target, mass rate limit violated
        ctrl = ProducerControl(BottomHolePressureTarget(50e5))
        new_ctrl, changed = JutulDarcy.check_well_limits(limits, limit_test_cond(50e5, -150.0), ctrl)
        @test changed
        @test new_ctrl.target isa TotalMassRateTarget
        @test new_ctrl.target.value == -100.0
    end
    @testset "Injector" begin
        limits = (bhp = 200e5, mrat = 100.0)
        # Mass rate target, BHP limit violated
        ctrl = InjectorControl(TotalMassRateTarget(100.0), [1.0]; density = 1000.0)
        new_ctrl, changed = JutulDarcy.check_well_limits(limits, limit_test_cond(250e5, 100.0), ctrl)
        @test changed
        @test new_ctrl.target isa BottomHolePressureTarget
        @test new_ctrl.target.value == 200e5
        # Mass rate target, BHP limit satisfied
        new_ctrl, changed = JutulDarcy.check_well_limits(limits, limit_test_cond(150e5, 100.0), ctrl)
        @test !changed
        # BHP target, mass rate limit violated
        ctrl = InjectorControl(BottomHolePressureTarget(200e5), [1.0]; density = 1000.0)
        new_ctrl, changed = JutulDarcy.check_well_limits(limits, limit_test_cond(200e5, 150.0), ctrl)
        @test changed
        @test new_ctrl.target isa TotalMassRateTarget
        @test new_ctrl.target.value == 100.0
    end
end
