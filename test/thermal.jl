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

    sys = CompositeSystem(:Reservoir, flow = sys_f, thermal = sys_t)
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
