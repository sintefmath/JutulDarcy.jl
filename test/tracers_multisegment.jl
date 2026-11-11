using Jutul, JutulDarcy, Test

@testset "Multisegment tracer fluxes" begin
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
    tracers = [SinglePhaseTracer(system, AqueousPhase()), MultiPhaseTracer(system)]
    flux = JutulDarcy.Tracers.TracerFluxType(tracers)
    concentrations = [0.5 0.75; 0.2 0.4]
    mass = [3.0 4.0; 5.0 6.0]
    state0 = (TracerMasses = zeros(2, 2),)
    state = (TracerMasses = mass, TracerConcentrations = concentrations,
        Saturations = [0.2 0.8; 0.8 0.2],
        PhaseMassDensities = [1000.0 1000.0; 800.0 800.0], TotalMassFlux = [10.0])

    function residual(state, connectivity = reshape([1, 2], 2, 1))
        well = MultiSegmentWell(connectivity, [1, 2], [1, 2])
        flow = JutulDarcy.WellSegmentFlow(well)
        eq = ConservationLaw(flow, :TracerMasses, 2; flux)
        out = zeros(2, 2)
        for cell in 1:2
            Jutul.update_equation_in_entity!(view(out, :, cell), cell,
                state, state0, eq, (system = system,), 2.0)
        end
        return out
    end

    forward_flux = [10*0.5*200/840, 10*0.2]
    forward = residual(state)
    @test forward ≈ mass/2 + hcat(forward_flux, -forward_flux)
    @test vec(sum(forward; dims = 2)) ≈ vec(sum(mass/2; dims = 2))
    reverse_state = merge(state, (TotalMassFlux = [-10.0],))
    reverse_flux = [-10*0.75*800/960, -10*0.4]
    @test residual(reverse_state) ≈ mass/2 + hcat(reverse_flux, -reverse_flux)
    @test residual(reverse_state, reshape([2, 1], 2, 1)) ≈ forward
    dry = merge(state, (Saturations = [0.0 0.8; 1.0 0.2],))
    @test residual(dry)[1, :] ≈ mass[1, :]/2
    stagnant = merge(state, (TotalMassFlux = [0.0],))
    @test residual(stagnant) ≈ mass/2
end

@testset "Multisegment tracer simulation" begin
    grid = CartesianMesh((3, 3, 1), (30.0, 30.0, 10.0))
    domain = reservoir_domain(grid; permeability = 100*si_unit(:milli)*si_unit(:darcy),
        porosity = 0.2)
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase());
        reference_densities = [1000.0, 800.0])
    injector = setup_well(domain, [1, 4, 7]; name = :I, simple_well = false)
    producer = setup_well(domain, [3, 6, 9]; name = :P, simple_well = false)
    model = setup_reservoir_model(domain, system; wells = [injector, producer],
        tracers = [SinglePhaseTracer(system, AqueousPhase()), MultiPhaseTracer(system)])
    state0 = setup_reservoir_state(model; Pressure = 200*si_unit(:bar),
        Saturations = [0.8, 0.2])
    controls = Dict(
        :I => InjectorControl(TotalRateTarget(1e-3), [1.0, 0.0]; density = 1000.0,
            tracers = [3.125, 1.0]),
        :P => ProducerControl(BottomHolePressureTarget(190*si_unit(:bar)))
    )
    forces = setup_reservoir_forces(model; control = controls)
    case = JutulCase(model, [60.0, 3600.0], forces; state0)
    result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true)
    @test length(result.states) == 2
    for state in result.states
        c = state[:TracerConcentrations]
        @test all(isfinite, c)
        @test minimum(c) >= -1e-8
    end
    @test all(i -> maximum(result.states[end][:TracerConcentrations][i, :]) > 0, 1:2)
    for well in (:I, :P)
        c = result.result.states[end][well][:TracerConcentrations]
        @test size(c) == (2, 3)
        @test all(isfinite, c)
    end
    @test minimum(result.result.states[end][:I][:TracerConcentrations]) > 0
end
