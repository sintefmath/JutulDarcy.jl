using Jutul, JutulDarcy, Test
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
    sys_f = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    sys_t = ThermalSystem(nphases = 2)

    sys = reservoir_system(flow = sys_f, thermal = sys_t)
    D = discretized_domain_tpfv_flow(G)
    if use_blocks
        l = BlockMajorLayout()
    else
        l = EquationMajorLayout()
    end
    ctx = DefaultContext(matrix_layout = l)

    model = SimulationModel(D, sys, data_domain = G, context = ctx)
    push!(model.output_variables, :Temperature)
    kr = BrooksCoreyRelativePermeabilities(2, [2.0, 2.0])
    replace_variables!(model, RelativePermeabilities = Pair(:flow, kr))
    forces_f = nothing
    forces = setup_forces(model, flow = forces_f)

    parameters = setup_parameters(model,
                                PhaseViscosities = [1e-3, 1e-3],
                                RockThermalConductivities = 1e-2,
                                FluidThermalConductivities = 1e-2,
                                RockDensity = 1e3,
                                FluidHeatCapacity = 10000.0,
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
using Jutul, JutulDarcy, GLMakie
function solve_thermal_wells(;
        nx = 10,
        ny = nx,
        nz = 2,
        block_backend = false,
        thermal = true,
        simple_well = false,
        composite = thermal,
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
        sys_f = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
        nph = 1
        i_mix = [1.0]
        rhoS = [rhoWS]
        c = [1e-6/bar]
    else
        # Set up a two-phase immiscible system
        phases = (AqueousPhase(), VaporPhase())
        rhoS = [rhoWS, rhoGS]
        sys_f = ImmiscibleSystem(phases, reference_densities = rhoS)
        nph = 2
        i_mix = [0.0, 1.0]
        c = [1e-6/bar, 1e-5/bar]
    end
    sys_t = ThermalSystem(nphases = nph)
    if composite
        if thermal
            sys = reservoir_system(flow = sys_f, thermal = sys_t)
        else
            sys = reservoir_system(flow = sys_f)
        end
    else
        if thermal
            sys = sys_t
        else
            sys = sys_f
        end
    end
    wells = [I, P]
    model, parameters = setup_reservoir_model(res, sys, wells = wells, block_backend = block_backend)
    # Replace the density function with our custom version for wells and reservoir
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    if composite
        tmp = Pair(:flow, ρ)
    else
        tmp = ρ
    end
    replace_variables!(model, PhaseMassDensities = tmp)
    dt = repeat([30.0]*day, 12*5)
    rate_target = TotalRateTarget(sum(parameters[:Reservoir][:FluidVolume])/sum(dt))
    bhp_target = BottomHolePressureTarget(50*bar)

    ictrl = InjectorControl(rate_target, i_mix, density = rhoGS, temperature = 300.0)

    pctrl = ProducerControl(bhp_target)

    controls = Dict(:Injector => ictrl,
                    :Producer => pctrl)
    # Wrap forces and initialize the state
    il = 3
    max_cuts = 0
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
    result = simulate_reservoir(state0, model, dt, forces = forces, parameters = parameters)
    @test length(result.states) == length(dt)
    return (result.states, result.result.reports, missing)
end
##

@testset "thermal wells" begin
    @testset "basic composite system" begin
        # Check that composite system gives same result before adding thermal
        states, reports, sim = solve_thermal_wells(nx = 3, thermal = false, composite = false);
        states_c, reports_c, sim_c = solve_thermal_wells(nx = 3, thermal = false, composite = true);
        using Test
        rstate = states[end]
        rstate_c = states_c[end]
        @test norm(rstate[:Pressure]-rstate_c[:Pressure])/norm(rstate[:Pressure]) < 1e-8
        @test norm(rstate[:Saturations]-rstate_c[:Saturations])/norm(rstate[:Saturations]) < 1e-6
    end
    @testset "well types" begin
        solve_thermal_wells(nx = 3, nz = 1,
            thermal = true, composite = true, simple_well = false, block_backend = false);
        solve_thermal_wells(nx = 3, nz = 1,
            thermal = true, composite = true, simple_well = true, block_backend = false);
    end
end

