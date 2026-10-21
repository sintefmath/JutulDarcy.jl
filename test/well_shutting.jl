using Jutul, JutulDarcy, Test

function update_shutting_test_rate(ctrl, q, dq; can_shut = true)
    cfg = JutulDarcy.WellGroupConfiguration([:Well])
    cfg.operating_controls[:Well] = ctrl
    cfg.requested_controls[:Well] = ctrl
    state = (TotalSurfaceMassRate = [q], WellGroupConfiguration = cfg)
    model = (domain = WellGroup([:Well]; can_shut_wells = can_shut),)
    Jutul.update_primary_variable!(state, TotalSurfaceMassRate(),
        :TotalSurfaceMassRate, model, [dq], 1.0)
    return state, cfg, model
end

@testset "Zero-rate runtime well closure" begin
    producer = ProducerControl(BottomHolePressureTarget(180e5); factor = 0.8)
    injector = InjectorControl(BottomHolePressureTarget(220e5), [1.0, 0.0];
        density = 1000.0, temperature = 320.0, factor = 0.7)
    floor_rate = JutulDarcy.MIN_INITIAL_WELL_RATE
    for (ctrl, sign) in ((producer, -1.0), (injector, 1.0))
        @testset "$(typeof(ctrl))" begin
            # The shutdown decision must depend on primal rates at the floor,
            # regardless of the sign of the AD partial.
            for tangent in (-1.0, 1.0), can_shut in (false, true)
                q = tangent*Jutul.get_ad_entity_scalar(
                    tangent*sign*floor_rate, 1, 1; tag = JutulDarcy.Wells())
                state, cfg, model = update_shutting_test_rate(ctrl, q,
                    -2sign*floor_rate; can_shut)
                operating = cfg.operating_controls[:Well]
                @test cfg.requested_controls[:Well] === ctrl
                @test JutulDarcy.ForwardDiff.partials(only(state.TotalSurfaceMassRate)) ==
                    JutulDarcy.ForwardDiff.partials(q)
                if !can_shut
                    @test operating === ctrl
                    @test Jutul.value(only(state.TotalSurfaceMassRate)) == sign*floor_rate
                    continue
                end
                @test JutulDarcy.is_zero_rate_control(operating)
                @test operating isa (ctrl isa ProducerControl ? ProducerControl : InjectorControl)
                @test operating.target isa TotalRateTarget
                @test iszero(operating.target.value)
                @test operating.factor == ctrl.factor
                @test iszero(Jutul.value(only(state.TotalSurfaceMassRate)))
                if ctrl isa InjectorControl
                    @test operating.injection_mixture == ctrl.injection_mixture
                    @test operating.mixture_density == ctrl.mixture_density
                    @test operating.temperature == ctrl.temperature
                end

                # A subsequent Newton update must keep the rate at zero,
                # rather than restoring the nonzero startup floor.
                Jutul.update_primary_variable!(state, TotalSurfaceMassRate(),
                    :TotalSurfaceMassRate, model, [sign*1e-3], 1.0)
                @test iszero(Jutul.value(only(state.TotalSurfaceMassRate)))
                @test iszero(Jutul.value(JutulDarcy.valid_surface_rate_for_control(q, operating)))

                # An infeasible BHP is expected after runtime closure. It must
                # neither reactivate the BHP target nor enter complementarity.
                cfg.limits[:Well] = (bhp = ctrl.target.value,)
                bhp = ctrl.target.value + sign*10e5
                cond = JutulDarcy.FacilityVariablesForWell(1, :Well,
                    bhp, 0.0, 0.0, 0.0, 0.0)
                JutulDarcy.apply_well_limits!(cfg, model, state,
                    cfg.limits[:Well], operating, :Well, cond)
                @test cfg.operating_controls[:Well] === operating
                @test iszero(JutulDarcy.well_control_equation(
                    operating, cond, :Well, model, state))

                # The zero-rate equation remains linear with a nonzero
                # derivative, including at exactly zero surface rate.
                qw = Jutul.get_ad_entity_scalar(0.0, 1, 1; tag = JutulDarcy.Wells())
                ad_cond = JutulDarcy.FacilityVariablesForWell(1, :Well,
                    bhp, 0.0, qw, 0.0, 0.0)
                residual = JutulDarcy.well_control_equation(
                    operating, ad_cond, :Well, model, state)
                @test iszero(Jutul.value(residual))
                @test only(JutulDarcy.ForwardDiff.partials(residual)) ≈ 86400.0
            end
        end
    end
end
