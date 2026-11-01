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
            model = SimulationModel(domain, JutulDarcy.FacilitySystem(sys))
            controls = Dict{Symbol, Any}(w => control(target) for w in wells)
            controls[:D] = control(TotalRateTarget(sgn*7e-4))
            controls[:Off] = DisabledControl()
            cfg = JutulDarcy.WellGroupConfiguration(wells, controls)
            rates = sgn*repeat([1e-4, 2e-4, 4e-4], 1, length(wells))
            state = JutulStorage(Dict(
                :BottomHolePressure => [1e7, 1e7, 2e7, 1e7, 1e7],
                :TotalSurfaceMassRate => vec(sum(rates; dims = 1)),
                :SurfacePhaseRates => rates,
                :WellGroupConfiguration => cfg
            ))
            remaining = target.value - sum(rates[ph, 4] for ph in phases)
            fractions = (1/2, 1/2)
            if potentials
                fractions = (1/3, 2/3)
            end
            for (well, fraction) in zip((:B, :C), fractions)
                cond = JutulDarcy.FacilityVariablesForWell(model, state, well)
                @test JutulDarcy.get_control_target_value(controls[well], target, cond, well, model, state) ≈ fraction*remaining
                @test JutulDarcy.well_target_value(controls[well], target, cond, well, model, state) ≈ sum(rates[ph, cond.idx] for ph in phases)
            end
            @test JutulDarcy.as_limit(target) == NamedTuple{(symbol,)}((target.value,))
            @test JutulDarcy.translate_target_to_symbol(target) == symbol
            @test JutulDarcy.setup_injector_control(0.01, symbol, [1.0, 0.0, 0.0]).target isa Target
            @test JutulDarcy.setup_producer_control(0.01, symbol).target.value == -0.01
            @test_throws ArgumentError Target(Inf)
            @test_throws ArgumentError Target(NaN)
            @test Target(Float32(0.01)).value isa Float32
            @test occursin(string(Target), sprint(show, target))

            # AD partials must not trigger switching at the allocated rate.
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            allocated = JutulDarcy.get_control_target_value(controls[:B], target, cond, :B, model, state)
            boundary_rates = zeros(3)
            boundary_rates[first(phases)] = allocated
            boundary = JutulDarcy.FacilityVariablesForWell(3, :B, cond.bottom_hole_pressure,
                cond.total_mass_rate, boundary_rates...)
            for tangent in (-1.0, 1.0)
                ad_bhp = [tangent*Jutul.get_ad_entity_scalar(tangent*p, 1, 1; tag = Wells())
                    for p in state.BottomHolePressure]
                ad_state = JutulStorage(Dict(:BottomHolePressure => ad_bhp,
                    :TotalSurfaceMassRate => state.TotalSurfaceMassRate,
                    :SurfacePhaseRates => state.SurfacePhaseRates, :WellGroupConfiguration => cfg))
                @test ismissing(JutulDarcy.check_well_limit(symbol, target.value,
                    boundary, controls[:B], :B, model, ad_state))
            end

            # Leave group control at a BHP limit, then rejoin when the group rate binds.
            bhp_limit = 2.5e7
            if inject
                bhp_limit = 1.5e7
            end
            limits = merge(JutulDarcy.as_limit(target), (bhp = bhp_limit,))
            cfg.limits[:B] = limits
            # A group residual must not prevent switching to an individual limit.
            # Test both a nearly converged allocation and a large Newton residual.
            saved_rates = copy(state.SurfacePhaseRates[:, 3])
            for excess in (1e-7, 1.0)
                cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
                allocated = JutulDarcy.get_control_target_value(controls[:B], target, cond, :B, model, state)
                state.SurfacePhaseRates[:, 3] .= 0.0
                state.SurfacePhaseRates[first(phases), 3] = (1 + excess)*allocated
                cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
                @test !ismissing(JutulDarcy.check_well_limit(symbol, target.value,
                    cond, controls[:B], :B, model, state))
                switched, changed = JutulDarcy.check_well_limits(limits, cond, controls[:B], :B, model, state)
                @test changed
                @test switched.target == BottomHolePressureTarget(bhp_limit)
            end
            state.SurfacePhaseRates[:, 3] .= saved_rates
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            JutulDarcy.apply_well_limits!(cfg, model, state, limits, controls[:B], :B, cond)
            limited = cfg.operating_controls[:B]
            @test limited.target isa BottomHolePressureTarget
            @test limited.target.value == bhp_limit

            state.BottomHolePressure[3] = bhp_limit
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            allocated = JutulDarcy.get_control_target_value(control(target), target, cond, :B, model, state)
            state.SurfacePhaseRates[:, 3] .= 0.0
            state.SurfacePhaseRates[first(phases), 3] = 2allocated
            cond = JutulDarcy.FacilityVariablesForWell(model, state, :B)
            JutulDarcy.apply_well_limits!(cfg, model, state, limits, limited, :B, cond)
            @test cfg.operating_controls[:B].target == target
            @test cfg.requested_controls[:B].target == target
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
    @test model[:Facility].domain.group_to_well_numbers[:Production] == [4, 2]
    state0 = setup_reservoir_state(model; Pressure = 200e5, Saturations = [0.4, 0.4, 0.2])
    rate = 1e-3
    volume_fractions = [0.4, 0.4, 0.2]
    density = sum(volume_fractions.*densities)
    mixture = volume_fractions.*densities/density
    injector = InjectorControl(Target(rate), mixture; density,
        phases = Tuple(enumerate(volume_fractions)))
    producer = ProducerControl(Target(-rate))
    controls = Dict(:I1 => injector, :I2 => injector, :P1 => producer, :P2 => producer)
    # One well in each group is rate limited; the other takes the remainder.
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
            JutulDarcy.individual_rate_target(Target(0.999rate)))
        next_controls[:P2] = ProducerControl(JutulDarcy.individual_rate_target(Target(-0.999rate)))
        next_forces = setup_reservoir_forces(model; control = next_controls, limits)
        forces = [forces, next_forces]
    end
    result = simulate_reservoir(state0, model, fill(0.1*si_unit(:day), 2);
        forces, info_level = -1, timesteps = :none, error_on_incomplete = true,
        nonlinear_tolerance = 1e-10, tolerances = (tol_control = 1e-10,))
    @test length(result.states) == 2
    for (step, state) in enumerate(result.result.states)
        facility = state[:Facility]
        rates = facility[:SurfacePhaseRates]
        cfg = facility[:WellGroupConfiguration]
        @test cfg.limits[:I1][symbol] == rate
        @test cfg.limits[:P1][symbol] == -rate
        @test sum(rates[ph, w] for ph in phases, w in (1, 3)) ≈ rate rtol = 1e-5
        @test sum(rates[ph, w] for ph in phases, w in (2, 4)) ≈ -rate rtol = 1e-5
        if switch_back && step == 2
            @test cfg.operating_controls[:I1].target isa Target
            @test cfg.operating_controls[:P1].target isa Target
            @test sum(rates[ph, 1] for ph in phases) ≈ 0.001rate rtol = 1e-5
            @test sum(rates[ph, 2] for ph in phases) ≈ -0.001rate rtol = 1e-5
        elseif bhp_limits
            @test facility[:BottomHolePressure][1] ≈ limits[:I1].bhp rtol = 1e-6
            @test facility[:BottomHolePressure][2] ≈ limits[:P1].bhp rtol = 1e-6
            @test sum(rates[:, 1]) < rate/2
            @test -sum(rates[:, 2]) < rate/2
        else
            @test sum(rates[:, 1]) ≈ 2e-4 rtol = 1e-5
            @test sum(rates[:, 2]) ≈ -2e-4 rtol = 1e-5
        end
        @test all(rates[:, [1, 3]] .>= 0.0)
        @test all(rates[:, [2, 4]] .<= 0.0)
    end
end

@testset "Injection and production groups with limits" begin
    for (Target, symbol, phases) in group_rate_cases, potentials in (false, true)
        @testset "$symbol, potentials=$potentials" begin
            well_group_simulation(Target, phases, potentials)
        end
    end
    for potentials in (false, true)
        @testset "BHP limits, potentials=$potentials" begin
            well_group_simulation(GroupTotalRateTarget, (1, 2, 3), potentials; bhp_limits = true)
        end
        @testset "Return to group control, potentials=$potentials" begin
            well_group_simulation(GroupTotalRateTarget, (1, 2, 3), potentials; switch_back = true)
        end
    end
end
