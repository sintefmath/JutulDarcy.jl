using Jutul, JutulDarcy, Test

function setup_two_cell_well_shutting(is_injector; operable = true)
    bar = si_unit(:bar)
    mesh = CartesianMesh((2, 1, 1), (20.0, 10.0, 10.0))
    domain = reservoir_domain(mesh, permeability = si_unit(:darcy), porosity = 0.2)
    well = setup_well(domain, [1], name = :Well, simple_well = true)
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, sys, wells = [well], block_backend = false)
    state0 = setup_reservoir_state(model, Pressure = 200bar, Saturations = [0.8, 0.2])
    # Cell two is connected to a constant-pressure aquifer. Reversing the BHP
    # limit makes the well inoperable regardless of its requested rate.
    bhp = ((is_injector == operable) ? 250.0 : 150.0)*bar
    rate = is_injector ? 1e-4 : -1.0
    if is_injector
        control = InjectorControl(BottomHolePressureTarget(bhp), [1.0, 0.0],
            density = 1000.0, temperature = 320.0, enthalpy = 1e5,
            tracers = [0.25], factor = 0.8)
    else
        control = ProducerControl(TotalRateTarget(rate), factor = 0.8)
    end
    forces = setup_reservoir_forces(model,
        control = Dict(:Well => control),
        limits = Dict(:Well => (bhp = bhp, rate = rate)),
        set_default_limits = false, bc = FlowBoundaryCondition(2, 200bar))
    sim, config = setup_reservoir_simulator(model, state0,
        info_level = -1, linear_solver = nothing, failure_cuts_timestep = false)
    # Prepare explicitly so the tests can prescribe the Newton update that
    # reaches the shutting condition, then continue the coupled flow solve.
    Jutul.update_secondary_variables!(sim.storage, sim.model)
    Jutul.update_before_step!(sim, si_unit(:day), forces)
    return sim, config, forces
end

function update_towards_wrong_well_direction!(sim, is_injector)
    state = sim.storage.Facility.state
    dx = [(is_injector ? -2.0 : 2.0)*JutulDarcy.MIN_INITIAL_WELL_RATE]
    Jutul.update_primary_variable!(state, JutulDarcy.TotalSurfaceMassRate(),
        :TotalSurfaceMassRate, sim.model[:Facility], dx, 1.0)
end

function finish_two_cell_well_step!(sim, config, forces)
    done, report = Jutul.solve_ministep(sim, si_unit(:day), forces, 15, config, prepare = false)
    @test done
    @test report[:success]
    @test all(isfinite, Jutul.value.(sim.storage.Reservoir.state.Pressure))
    @test all(isfinite, Jutul.value.(sim.storage.Facility.state.BottomHolePressure))
end

@testset "Two-cell well shutting" begin
    for is_injector in (false, true)
        @testset "$(is_injector ? "Injector" : "Producer")" begin
            @testset "Try alternative control before shutting" begin
                sim, config, forces = setup_two_cell_well_shutting(is_injector)
                state = sim.storage.Facility.state
                cfg = state.WellGroupConfiguration
                requested = cfg.requested_controls[:Well]
                initial_rate = copy(state.TotalSurfaceMassRate)
                # Use the actual AD state at the rate floor: derivatives must
                # not affect the decision to try the alternative control.
                update_towards_wrong_well_direction!(sim, is_injector)
                control = cfg.operating_controls[:Well]
                if is_injector
                    @test control.target isa TotalRateTarget
                    @test control.target.value == cfg.limits[:Well].rate
                    for field in (:injection_mixture, :mixture_density, :phases,
                            :temperature, :enthalpy, :tracers)
                        @test isequal(getproperty(control, field), getproperty(requested, field))
                    end
                else
                    @test control.target isa BottomHolePressureTarget
                    @test control.target.value == cfg.limits[:Well].bhp
                end
                @test control.factor == requested.factor
                @test cfg.requested_controls[:Well] === requested
                @test state.TotalSurfaceMassRate == initial_rate
                finish_two_cell_well_step!(sim, config, forces)
                @test !(cfg.operating_controls[:Well] isa DisabledControl)
                rate = only(Jutul.value.(state.TotalSurfaceMassRate))
                if is_injector
                    @test rate ≈ requested.mixture_density*cfg.limits[:Well].rate
                else
                    @test rate < 0.0
                    @test only(Jutul.value.(state.BottomHolePressure)) ≈ cfg.limits[:Well].bhp
                end
            end

            @testset "Shut if the alternative also fails" begin
                sim, config, forces = setup_two_cell_well_shutting(is_injector)
                state = sim.storage.Facility.state
                cfg = state.WellGroupConfiguration
                update_towards_wrong_well_direction!(sim, is_injector)
                update_towards_wrong_well_direction!(sim, is_injector)
                @test cfg.operating_controls[:Well] isa DisabledControl
                @test iszero(only(Jutul.value.(state.TotalSurfaceMassRate)))
                finish_two_cell_well_step!(sim, config, forces)
                @test cfg.operating_controls[:Well] isa DisabledControl
                @test all(iszero, Jutul.value.(state.SurfacePhaseRates))

                # Shutting and the record of the attempted alternative are
                # temporary: the well gets a fresh attempt in the next step.
                Jutul.update_before_step!(sim, si_unit(:day), forces)
                update_towards_wrong_well_direction!(sim, is_injector)
                @test !(cfg.operating_controls[:Well] isa DisabledControl)
                finish_two_cell_well_step!(sim, config, forces)
            end

            @testset "Infeasible BHP shuts without cycling" begin
                sim, config, forces = setup_two_cell_well_shutting(is_injector, operable = false)
                state = sim.storage.Facility.state
                targets = Any[]
                config[:post_iteration_hook] = (converged, report, storage, model, dt, f, c, it) -> begin
                    push!(targets, storage.Facility.state.WellGroupConfiguration.operating_controls[:Well].target)
                    return converged
                end
                finish_two_cell_well_step!(sim, config, forces)
                @test state.WellGroupConfiguration.operating_controls[:Well] isa DisabledControl
                @test iszero(only(Jutul.value.(state.TotalSurfaceMassRate)))
                @test all(iszero, Jutul.value.(state.SurfacePhaseRates))
                if is_injector
                    @test any(t -> t isa TotalRateTarget, targets)
                end
            end

            @testset "No alternative limit" begin
                for limits in (nothing, NamedTuple())
                    sim, config, forces = setup_two_cell_well_shutting(is_injector)
                    state = sim.storage.Facility.state
                    state.WellGroupConfiguration.limits[:Well] = limits
                    update_towards_wrong_well_direction!(sim, is_injector)
                    @test state.WellGroupConfiguration.operating_controls[:Well] isa DisabledControl
                    finish_two_cell_well_step!(sim, config, forces)
                    @test iszero(only(Jutul.value.(state.TotalSurfaceMassRate)))
                end
            end
        end
    end
end
