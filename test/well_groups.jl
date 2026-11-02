using Jutul, JutulDarcy, Test

const group_rate_cases = (
    (GroupTotalRateTarget, :group_rate, (1, 2, 3)),
    (GroupSurfaceLiquidRateTarget, :group_lrat, (1, 2)),
    (GroupSurfaceOilRateTarget, :group_orat, (2,)),
    (GroupSurfaceGasRateTarget, :group_grat, (3,)),
    (GroupSurfaceWaterRateTarget, :group_wrat, (1,))
)

@testset "Well group construction" begin
    group = WellGroup([:P, :I2, :I1]; groups = Dict(:Injection => [:I1, :I2]))
    @test group.group_to_well_numbers[:Injection] == [3, 2]
    @test group.use_potentials == Dict(:Injection => true)
    @test !group.lag_potentials
    @test WellGroup([:W]; lag_potentials = true).lag_potentials
    @test_throws ArgumentError WellGroup([:W]; groups = Dict(:G => [:Unknown]))
    @test_throws ArgumentError WellGroup([:W]; groups = Dict(:G => [:W, :W]))
    @test_throws ArgumentError WellGroup([:W]; groups = Dict(:G1 => [:W], :G2 => [:W]))
    @test_throws ArgumentError WellGroup([:W]; groups = Dict(:G => Symbol[]))
    @test_throws ArgumentError WellGroup([:W]; use_potentials = Dict(:Unknown => false))
end

@testset "Group phase rates and allocation" begin
    for (Target, symbol, phases) in group_rate_cases, inject in (false, true), potentials in (false, true)
        @testset "$symbol, injection=$inject, potentials=$potentials" begin
            sgn = -1.0
            if inject
                sgn = 1.0
            end
            target = Target(sgn*0.01)
            function control(t)
                if inject
                    return InjectorControl(t, [1/3, 1/3, 1/3]; density = 700.0)
                else
                    return ProducerControl(t)
                end
            end
            wells = [:Outside, :C, :B, :D, :Off]
            domain = WellGroup(wells; groups = Dict(:G => [:B, :D, :C, :Off]),
                use_potentials = Dict(:G => potentials))
            sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase()))
            densities = JutulDarcy.reference_densities(sys)
            model = SimulationModel(domain, JutulDarcy.FacilitySystem(sys))
            controls = Dict{Symbol, Any}(w => control(target) for w in wells)
            controls[:D] = control(TotalRateTarget(sgn*7e-4))
            controls[:Off] = DisabledControl()
            cfg = JutulDarcy.WellGroupConfiguration(wells, controls)
            rates = sgn*repeat([1e-4, 2e-4, 4e-4], 1, length(wells))
            rates[:, 2] .*= 2.0
            rates[:, [1, 5]] .*= 10.0
            state = JutulStorage(Dict(
                :BottomHolePressure => [1e7, 2e7, 2e7, 1e7, 1e7],
                :TotalSurfaceMassRate => vec(sum(densities.*rates; dims = 1)),
                :SurfacePhaseRates => rates,
                :WellGroupConfiguration => cfg
            ))
            remaining = target.value - sum(rates[ph, 4] for ph in phases)
            allocations = Float64[]
            for well in (:B, :C)
                cond = JutulDarcy.FacilityVariablesForWell(model, state, well)
                push!(allocations, JutulDarcy.get_control_target_value(controls[well], target, cond, well, model, state))
            end
            # Group members share the remainder after individual limits. Neither
            # an unrelated well nor a disabled group member contributes.
            @test sum(allocations) ≈ remaining
            @test all(sgn.*allocations .> 0.0)
            if potentials
                # B and C have equal BHP, while C has twice B's phase rates.
                @test allocations ≈ remaining.*[1/3, 2/3]
            else
                @test allocations ≈ remaining.*[1/2, 1/2]
            end
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            @test JutulDarcy.well_target_value(controls[:B], target, cond, :B, model, state) ≈ sum(rates[ph, 3] for ph in phases)

            # Exceed the entire group rate, so that the rate limit is violated
            # even when the potential weights depend on the current rates.
            state.SurfacePhaseRates[:, 3] .= 0.0
            state.SurfacePhaseRates[first(phases), 3] = 2target.value
            state.TotalSurfaceMassRate[3] = sum(densities.*state.SurfacePhaseRates[:, 3])
            bhp_limit = 2.5e7
            if inject
                bhp_limit = 1.5e7
            end
            limits = merge(JutulDarcy.as_limit(target), (bhp = bhp_limit,))
            cfg.limits[:B] = limits
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            JutulDarcy.apply_well_limits!(cfg, model, state, limits, controls[:B], :B, cond)
            limited = cfg.operating_controls[:B]
            @test limited.target == BottomHolePressureTarget(bhp_limit)

            # At the BHP limit, excess production/injection must switch back to
            # group control. This also checks that the group target is a limit.
            state.BottomHolePressure[3] = bhp_limit
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            JutulDarcy.apply_well_limits!(cfg, model, state, limits, limited, :B, cond)
            @test cfg.operating_controls[:B].target == target
            @test cfg.requested_controls[:B].target == target
        end
    end
end

@testset "Group remainder at other phase constraints" begin
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase()))
    wells = [:P1, :P2]
    domain = WellGroup(wells; groups = Dict(:G => wells), use_potentials = Dict(:G => false))
    model = SimulationModel(domain, JutulDarcy.FacilitySystem(sys))
    controls = Dict(:P1 => ProducerControl(GroupSurfaceOilRateTarget(-0.01)),
        :P2 => ProducerControl(GroupSurfaceGasRateTarget(-0.02)))
    rates = -[0.001 0.002; 0.003 0.004; 0.005 0.006]
    cfg = JutulDarcy.WellGroupConfiguration(wells, controls)
    state = JutulStorage(Dict(:BottomHolePressure => fill(2e7, 2),
        :SurfacePhaseRates => rates,
        :TotalSurfaceMassRate => vec(sum(JutulDarcy.reference_densities(sys).*rates; dims = 1)),
        :WellGroupConfiguration => cfg))
    ctrl = controls[:P1]
    cond = JutulDarcy.FacilityVariablesForWell(model, state, :P1)
    allocated() = JutulDarcy.get_control_target_value(ctrl, ctrl.target, cond, :P1, model, state)
    @test allocated() ≈ -0.006
    state.SurfacePhaseRates[2, 2] = -0.02
    @test allocated() == 0.0
end

function well_group_test_case(Target = GroupTotalRateTarget; potentials = false, equal_densities = false)
    mesh = CartesianMesh((6, 4, 1), (600.0, 400.0, 30.0))
    permeability = [0.2e-12*(1 + i/24) for i in 1:24]
    domain = reservoir_domain(mesh; permeability, porosity = 0.2)
    # Interleave injectors/producers and reverse group order to exercise well lookup.
    wells = [setup_vertical_well(domain, i, j; name, simple_well = true)
        for (name, i, j) in ((:I1, 1, 1), (:P1, 6, 1), (:I2, 1, 4), (:P2, 6, 4))]
    groups = Dict(:Injection => [:I2, :I1], :Production => [:P2, :P1])
    densities = equal_densities ? fill(1000.0, 3) : [1000.0, 800.0, 100.0]
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase());
        reference_densities = densities)
    model = setup_reservoir_model(domain, sys; wells, well_groups = groups,
        group_use_potentials = Dict(:Injection => potentials, :Production => potentials))
    volume_fractions = [0.4, 0.4, 0.2]
    state0 = setup_reservoir_state(model; Pressure = 200e5, Saturations = volume_fractions)
    rate = 1e-3
    density = sum(volume_fractions.*densities)
    mixture = volume_fractions.*densities/density
    injector = InjectorControl(Target(rate), mixture; density,
        phases = Tuple(enumerate(volume_fractions)))
    producer = ProducerControl(Target(-rate))
    controls = Dict{Symbol, Any}(:I1 => injector, :I2 => injector, :P1 => producer, :P2 => producer)
    return (; model, state0, controls, rate)
end

function simulate_well_group_case(case, forces)
    nstep = forces isa AbstractVector ? length(forces) : 2
    result = simulate_reservoir(case.state0, case.model, fill(0.1*si_unit(:day), nstep);
        forces, info_level = -1, timesteps = :none, error_on_incomplete = true)
    @test length(result.result.states) == nstep
    return [state[:Facility] for state in result.result.states]
end

well_group_phase_rate(facility, phases, well_index) = sum(facility[:SurfacePhaseRates][ph, well_index] for ph in phases)

function test_group_totals(facility, phases, rate)
    # The absolute tolerance is 0.1% of the group target, also for small shares.
    rate_tol = 1e-3*rate
    @test sum(well_group_phase_rate(facility, phases, w) for w in (1, 3)) ≈ rate atol = rate_tol
    @test sum(well_group_phase_rate(facility, phases, w) for w in (2, 4)) ≈ -rate atol = rate_tol
    rates = facility[:SurfacePhaseRates]
    @test all(rates[:, [1, 3]] .>= 0.0)
    @test all(rates[:, [2, 4]] .<= 0.0)
end

@testset "Group control integration" begin
    @testset "Equal allocation and individual rate limits" begin
        for (Target, symbol, phases) in group_rate_cases
            @testset "$symbol" begin
                case = well_group_test_case(Target)
                (; model, controls, rate) = case
                unlimited = setup_reservoir_forces(model; control = controls)
                cap = 0.2rate
                limits = Dict(:I1 => (rate = cap, bhp = 300e5), :P1 => (rate = -cap, bhp = 100e5))
                limited = setup_reservoir_forces(model; control = controls, limits)
                @test limited[:Facility].limits[:I1][symbol] == rate
                @test limited[:Facility].limits[:P1][symbol] == -rate
                states = simulate_well_group_case(case, [unlimited, limited])
                for (step, facility) in enumerate(states)
                    test_group_totals(facility, phases, rate)
                    cfg = facility[:WellGroupConfiguration]
                    for (w, ctrl) in controls
                        @test cfg.requested_controls[w] == ctrl
                    end
                    if step == 1
                        # Equal weighting applies to the targeted phases, not total flow.
                        for (w, sgn) in ((1, 1), (3, 1), (2, -1), (4, -1))
                            @test well_group_phase_rate(facility, phases, w) ≈ sgn*rate/2 atol = 1e-3*rate
                        end
                    else
                        # Both caps must actually bind; the other wells take the remainder.
                        @test cfg.operating_controls[:I1].target == TotalRateTarget(cap)
                        @test cfg.operating_controls[:P1].target == TotalRateTarget(-cap)
                        @test cfg.operating_controls[:I2].target == Target(rate)
                        @test cfg.operating_controls[:P2].target == Target(-rate)
                        @test sum(facility[:SurfacePhaseRates][:, 1]) ≈ cap atol = 1e-3*rate
                        @test sum(facility[:SurfacePhaseRates][:, 2]) ≈ -cap atol = 1e-3*rate
                    end
                end
            end
        end
    end

    @testset "Potential weighting" begin
        case = well_group_test_case(GroupSurfaceWaterRateTarget; potentials = true)
        (; model, controls, rate) = case
        unlimited = setup_reservoir_forces(model; control = controls)
        cap = 0.2rate
        limits = Dict(:I1 => (rate = cap,), :P1 => (rate = -cap,))
        limited = setup_reservoir_forces(model; control = controls, limits)
        states = simulate_well_group_case(case, [unlimited, limited])
        for (step, facility) in enumerate(states)
            test_group_totals(facility, (1,), rate)
            if step == 1
                # Different well capacities must give different shares in both groups.
                @test abs(well_group_phase_rate(facility, (1,), 1) - rate/2) > 1e-3*rate
                @test abs(well_group_phase_rate(facility, (1,), 2) + rate/2) > 1e-3*rate
                @test sum(facility[:SurfacePhaseRates][:, 1]) > cap
                @test -sum(facility[:SurfacePhaseRates][:, 2]) > cap
            else
                # Potential weights change with the solution, so a well can end
                # up below its cap while the other member takes the remainder.
                @test sum(facility[:SurfacePhaseRates][:, 1]) <= cap + 1e-3*rate
                @test -sum(facility[:SurfacePhaseRates][:, 2]) <= cap + 1e-3*rate
            end
        end
    end

    @testset "BHP limits and return to group control" begin
        case = well_group_test_case(; equal_densities = true)
        (; model, controls, rate) = case
        limits = Dict(:I1 => (bhp = 200.5e5,), :P1 => (bhp = 199.5e5,))
        limited = setup_reservoir_forces(model; control = controls, limits)
        # Keep I1/P1's requested controls unchanged. Once I2/P2 take most of the
        # group rate, I1/P1 must automatically leave BHP control for the remainder.
        next_controls = copy(controls)
        next_controls[:I2] = JutulDarcy.replace_target(controls[:I2], TotalRateTarget(0.98rate))
        next_controls[:P2] = ProducerControl(TotalRateTarget(-0.98rate))
        next_forces = setup_reservoir_forces(model; control = next_controls, limits)
        states = simulate_well_group_case(case, [limited, next_forces, next_forces])
        for (step, facility) in enumerate(states)
            test_group_totals(facility, (1, 2, 3), rate)
            cfg = facility[:WellGroupConfiguration]
            @test cfg.requested_controls[:I1] == controls[:I1]
            @test cfg.requested_controls[:P1] == controls[:P1]
            if step == 1
                @test cfg.operating_controls[:I1].target == BottomHolePressureTarget(limits[:I1].bhp)
                @test cfg.operating_controls[:P1].target == BottomHolePressureTarget(limits[:P1].bhp)
                @test facility[:BottomHolePressure][1] ≈ limits[:I1].bhp atol = 0.01*si_unit(:bar)
                @test facility[:BottomHolePressure][2] ≈ limits[:P1].bhp atol = 0.01*si_unit(:bar)
                @test well_group_phase_rate(facility, (1, 2, 3), 1) < rate/2
                @test -well_group_phase_rate(facility, (1, 2, 3), 2) < rate/2
            else
                @test cfg.operating_controls[:I1].target == controls[:I1].target
                @test cfg.operating_controls[:P1].target == controls[:P1].target
                @test well_group_phase_rate(facility, (1, 2, 3), 1) ≈ 0.02rate atol = 1e-3*rate
                @test well_group_phase_rate(facility, (1, 2, 3), 2) ≈ -0.02rate atol = 1e-3*rate
                @test facility[:BottomHolePressure][1] < limits[:I1].bhp
                @test facility[:BottomHolePressure][2] > limits[:P1].bhp
            end
        end
    end

    @testset "Infeasible group target respects BHP limits" begin
        case = well_group_test_case(; equal_densities = true)
        (; model, controls, rate) = case
        limits = Dict(:I1 => (bhp = 200.5e5,), :I2 => (bhp = 200.5e5,),
            :P1 => (bhp = 199.5e5,), :P2 => (bhp = 199.5e5,))
        forces = setup_reservoir_forces(model; control = controls, limits)
        for facility in simulate_well_group_case(case, forces)
            cfg = facility[:WellGroupConfiguration]
            for (w, ix) in ((:I1, 1), (:P1, 2), (:I2, 3), (:P2, 4))
                @test cfg.operating_controls[w].target == BottomHolePressureTarget(limits[w].bhp)
                @test facility[:BottomHolePressure][ix] ≈ limits[w].bhp atol = 0.01*si_unit(:bar)
            end
            rates = facility[:SurfacePhaseRates]
            @test 0.0 < sum(rates[:, [1, 3]]) < 0.9rate
            @test 0.0 < -sum(rates[:, [2, 4]]) < 0.9rate
        end
    end

    @testset "Disabled members and reopening" begin
        case = well_group_test_case()
        (; model, controls, rate) = case
        disabled_controls = copy(controls)
        disabled_controls[:I1] = DisabledControl()
        disabled_controls[:P1] = DisabledControl()
        disabled = setup_reservoir_forces(model; control = disabled_controls)
        enabled = setup_reservoir_forces(model; control = controls)
        states = simulate_well_group_case(case, [disabled, enabled])
        for (step, facility) in enumerate(states)
            test_group_totals(facility, (1, 2, 3), rate)
            cfg = facility[:WellGroupConfiguration]
            if step == 1
                @test cfg.operating_controls[:I1] isa DisabledControl
                @test cfg.operating_controls[:P1] isa DisabledControl
                @test all(iszero, facility[:SurfacePhaseRates][:, [1, 2]])
            else
                @test cfg.operating_controls[:I1].target == controls[:I1].target
                @test cfg.operating_controls[:P1].target == controls[:P1].target
                @test well_group_phase_rate(facility, (1, 2, 3), 1) ≈ rate/2 atol = 1e-3*rate
                @test well_group_phase_rate(facility, (1, 2, 3), 2) ≈ -rate/2 atol = 1e-3*rate
            end
        end
    end
end
