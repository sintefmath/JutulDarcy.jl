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
    @test_throws ArgumentError WellGroup([:W]; groups = Dict(:G => [:Unknown]))
    @test_throws ArgumentError WellGroup([:W]; groups = Dict(:G => [:W, :W]))
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
                @test !isapprox(allocations[1], allocations[2])
            else
                @test allocations[1] ≈ allocations[2]
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
            @test limited.target isa BottomHolePressureTarget

            # At the BHP limit, excess production/injection must switch back to
            # group control. This also checks that the group target is a limit.
            state.BottomHolePressure[3] = bhp_limit
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            JutulDarcy.apply_well_limits!(cfg, model, state, limits, limited, :B, cond)
            @test cfg.operating_controls[:B].target == target
        end
    end
end

function well_group_simulation(Target, phases, potentials; bhp_limits = false, switch_back = false)
    mesh = CartesianMesh((6, 4, 1), (600.0, 400.0, 30.0))
    permeability = [0.2e-12*(1 + i/24) for i in 1:24]
    domain = reservoir_domain(mesh; permeability, porosity = 0.2)
    wells = [setup_vertical_well(domain, i, j; name, simple_well = true)
        for (name, i, j) in ((:I1, 1, 1), (:P1, 6, 1), (:I2, 1, 4), (:P2, 6, 4))]
    groups = Dict(:Injection => [:I2, :I1], :Production => [:P2, :P1])
    densities = [1000.0, 800.0, 100.0]
    if bhp_limits
        densities .= 1000.0
    end
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase());
        reference_densities = densities)
    model = setup_reservoir_model(domain, sys; wells, well_groups = groups,
        group_use_potentials = Dict(:Injection => potentials, :Production => potentials))
    state0 = setup_reservoir_state(model; Pressure = 200e5, Saturations = [0.4, 0.4, 0.2])
    rate = 1e-3
    volume_fractions = [0.4, 0.4, 0.2]
    density = sum(volume_fractions.*densities)
    mixture = volume_fractions.*densities/density
    injector = InjectorControl(Target(rate), mixture; density,
        phases = Tuple(enumerate(volume_fractions)))
    producer = ProducerControl(Target(-rate))
    controls = Dict(:I1 => injector, :I2 => injector, :P1 => producer, :P2 => producer)
    # Cap one well in each group; the other can take the remainder.
    limits = Dict(:I1 => (rate = 2e-4, bhp = 300e5), :P1 => (rate = -2e-4, bhp = 100e5))
    if bhp_limits
        limits = Dict(:I1 => (bhp = 200.5e5,), :P1 => (bhp = 199.5e5,))
    end
    forces = setup_reservoir_forces(model; control = controls, limits)
    symbol = JutulDarcy.translate_target_to_symbol(Target(rate))
    @test forces[:Facility].limits[:I1][symbol] == rate
    @test forces[:Facility].limits[:P1][symbol] == -rate
    if switch_back
        next_controls = copy(controls)
        next_controls[:I2] = JutulDarcy.replace_target(injector,
            JutulDarcy.individual_rate_target(Target(0.9rate)))
        next_controls[:P2] = ProducerControl(JutulDarcy.individual_rate_target(Target(-0.9rate)))
        next_forces = setup_reservoir_forces(model; control = next_controls, limits)
        forces = [forces, next_forces]
    end
    result = simulate_reservoir(state0, model, fill(0.1*si_unit(:day), 2);
        forces, info_level = -1, timesteps = :none, error_on_incomplete = true)
    @test length(result.states) == 2
    for (step, state) in enumerate(result.result.states)
        facility = state[:Facility]
        rates = facility[:SurfacePhaseRates]
        cfg = facility[:WellGroupConfiguration]
        # Rates are checked to 0.1% of the group target, including small remainders.
        rate_tol = 1e-3*rate
        @test sum(rates[ph, w] for ph in phases, w in (1, 3)) ≈ rate atol = rate_tol
        @test sum(rates[ph, w] for ph in phases, w in (2, 4)) ≈ -rate atol = rate_tol
        if switch_back && step == 2
            @test cfg.operating_controls[:I1].target isa Target
            @test cfg.operating_controls[:P1].target isa Target
            @test sum(rates[ph, 1] for ph in phases) ≈ 0.1rate atol = rate_tol
            @test sum(rates[ph, 2] for ph in phases) ≈ -0.1rate atol = rate_tol
        elseif bhp_limits
            @test facility[:BottomHolePressure][1] <= limits[:I1].bhp + 0.01*si_unit(:bar)
            @test facility[:BottomHolePressure][2] >= limits[:P1].bhp - 0.01*si_unit(:bar)
            if !potentials
                @test cfg.operating_controls[:I1].target isa BottomHolePressureTarget
                @test cfg.operating_controls[:P1].target isa BottomHolePressureTarget
            end
        else
            # A cap need not bind when potential weighting gives a well less.
            @test sum(rates[:, 1]) <= 2e-4 + rate_tol
            @test -sum(rates[:, 2]) <= 2e-4 + rate_tol
        end
        @test all(rates[:, [1, 3]] .>= 0.0)
        @test all(rates[:, [2, 4]] .<= 0.0)
    end
end

@testset "Injection and production groups with limits" begin
    for (Target, symbol, phases) in group_rate_cases
        @testset "$symbol" begin
            well_group_simulation(Target, phases, false)
        end
    end
    @testset "Potential weighting" begin
        well_group_simulation(GroupSurfaceWaterRateTarget, (1,), true)
    end
    @testset "BHP limits" begin
        well_group_simulation(GroupTotalRateTarget, (1, 2, 3), false; bhp_limits = true)
    end
    @testset "Return to group control" begin
        well_group_simulation(GroupTotalRateTarget, (1, 2, 3), false; switch_back = true)
    end
end
