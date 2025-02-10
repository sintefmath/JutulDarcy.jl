using Jutul, JutulDarcy, Test, LinearAlgebra
function solve_thermal(;
        nc = 10,
        time = 1000.0,
        nstep = 100,
        poro = 0.1,
        perm = 9.8692e-14,
        use_blocks = false
    )
    T = time
    tstep = repeat([T/nstep], nstep)
    G = get_1d_reservoir(nc, poro = poro, perm = perm)
    nc = number_of_cells(G)

    G[:porosity][1] *= 1000
    G[:porosity][end] *= 1000
    bar = 1e5
    p0 = repeat([1000*bar], nc)
    p0[1] = 2000*bar
    p0[end] = 500*bar

    s0 = zeros(2, nc)
    s0[2, :] .= 1.0
    s0[1, 1] = 1.0
    s0[1, 2] = 0.0

    # Define system and realize on grid
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    D = discretized_domain_tpfv_flow(G)
    if use_blocks
        l = BlockMajorLayout()
    else
        l = EquationMajorLayout()
    end
    ctx = DefaultContext(matrix_layout = l)

    model = SimulationModel(D, sys, data_domain = G, context = ctx)
    JutulDarcy.add_thermal_to_model!(model)
    push!(model.output_variables, :Temperature)
    kr = BrooksCoreyRelativePermeabilities(2, [2.0, 2.0])
    replace_variables!(model, RelativePermeabilities = kr)
    forces = setup_forces(model)

    parameters = setup_parameters(model,
                                PhaseViscosities = [1e-3, 1e-3],
                                RockThermalConductivities = 1e-2,
                                FluidThermalConductivities = 1e-2,
                                RockDensity = 1e3,
                                ComponentHeatCapacity = 10000.0,
                                RockHeatCapacity = 500.0)
    T0 = repeat([273.15], nc)
    T0[1] = 500.0
    state0 = setup_state(model, Pressure = p0, Saturations = s0, Temperature = T0)

    states, reports = simulate(state0, model, tstep,
        parameters = parameters, forces = forces, info_level = -1)
    pushfirst!(states, state0)
    return states, reports
end

using Test
@testset "simple_thermal" begin
    for use_blocks in [true, false]
        states, = solve_thermal(nc = 10, use_blocks = use_blocks);
        T = states[end][:Temperature]
        # Check that the first cell fullfills BC
        @test 490 < T[2] < 500
        # Check monotone temp profile
        @test all(x -> x < 0, diff(T))
    end
end
##
function solve_thermal_wells(;
        nx = 10,
        ny = nx,
        nz = 2,
        block_backend = false,
        thermal = true,
        simple_well = false,
        single_phase = false
    )
    day = 3600*24
    bar = 1e5
    g = CartesianMesh((nx, ny, nz), (1000.0, 1000.0, 100.0))
    Darcy = 9.869232667160130e-13
    K = repeat([0.1*Darcy], 1, number_of_cells(g))
    res = reservoir_domain(g, porosity = 0.1, permeability = K)
    # Vertical well in (1, 1, *), producer in (nx, ny, 1)
    P = setup_vertical_well(res, 1, 1, name = :Producer, simple_well = simple_well)
    I = setup_well(res, [(nx, ny, 1)], name = :Injector, simple_well = simple_well)
    rhoWS = 1000.0
    rhoGS = 700.0
    if single_phase
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
        nph = 1
        i_mix = [1.0]
        rhoS = [rhoWS]
        c = [1e-6/bar]
    else
        # Set up a two-phase immiscible system
        phases = (AqueousPhase(), VaporPhase())
        rhoS = [rhoWS, rhoGS]
        sys = ImmiscibleSystem(phases, reference_densities = rhoS)
        nph = 2
        i_mix = [0.0, 1.0]
        c = [1e-6/bar, 1e-5/bar]
    end
    wells = [I, P]
    model, parameters = setup_reservoir_model(res, sys,
        thermal = thermal,
        wells = wells,
        block_backend = block_backend
    )
    # Replace the density function with our custom version for wells and reservoir
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    replace_variables!(model, PhaseMassDensities = ρ)
    dt = repeat([30.0]*day, 12*5)
    rate_target = TotalRateTarget(sum(parameters[:Reservoir][:FluidVolume])/sum(dt))
    bhp_target = BottomHolePressureTarget(50*bar)

    ictrl = InjectorControl(rate_target, i_mix, density = rhoGS, temperature = 300.0)

    pctrl = ProducerControl(bhp_target)

    controls = Dict(:Injector => ictrl,
                    :Producer => pctrl)
    forces = setup_reservoir_forces(model, control = controls)
    if thermal
        state0 = setup_reservoir_state(model,
            Pressure = 150*bar,
            Saturations = [1.0, 0.0],
            Temperature = 300.0
        )
    else
        state0 = setup_reservoir_state(model,
            Pressure = 150*bar,
            Saturations = [1.0, 0.0]
        )
    end
    result = simulate_reservoir(state0, model, dt,
        forces = forces, parameters = parameters, info_level = -1)
    return (result.states, result.result.reports, dt, model)
end

@testset "thermal wells" begin
    @testset "well types and backends" begin
        for simple_well in [true, false]
            for block_backend in [false, true]
                if block_backend
                    bs = "block"
                else
                    bs = "scalar"
                end
                if simple_well
                    ws = "simple"
                else
                    ws = "ms"
                end
                @testset "$bs $ws well" begin
                    states, reports, dt, = solve_thermal_wells(nx = 3, nz = 1,
                        thermal = true,
                        simple_well = simple_well,
                        block_backend = block_backend
                    );
                    @test length(states) == length(dt)
                end
            end
        end
    end
end
