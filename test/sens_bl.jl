using Jutul, JutulDarcy, Test

function get_sens(model, state0, parameters, tstep, forces, G)
    sim = Simulator(model, state0 = state0, parameters = parameters)
    states, reports = simulate(sim, tstep, forces = forces, extra_timing = false, info_level = -1)

    grad_adj = Jutul.solve_adjoint_sensitivities(model, states, reports, G,
                    forces = forces, state0 = state0, parameters = parameters,
                    extra_timing = false, raw_output = false)

    grad_numeric = Dict{Symbol, Any}()
    for k in keys(grad_adj)
        grad_numeric[k] = Jutul.solve_numerical_sensitivities(model, states, reports, G, k,
                        forces = forces, state0 = state0, parameters = parameters,
                        epsilon = 1e-6)
    end
    return grad_adj, grad_numeric
end

function setup_bl(;nc = 100, time = 1.0, nstep = 100)
    T = time
    tstep = repeat([T/nstep], nstep)
    G = get_1d_reservoir(nc)
    nc = number_of_cells(G)
    timesteps = tstep*3600*24 # Convert time-steps from days to seconds

    bar = 1e5
    p0 = 100*bar
    # Define system and realize on grid
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    model = SimulationModel(G, sys)
    kr = BrooksCoreyRelPerm(sys, [2.0, 2.0], [0.2, 0.2])
    replace_variables!(model, RelativePermeabilities = kr)
    tot_time = sum(timesteps)
    irate = 500*sum(G.grid.pore_volumes)/tot_time
    src  = [SourceTerm(1, irate, fractional_flow = [0.8, 0.2]), 
            SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
    forces = setup_forces(model, sources = src)

    parameters = setup_parameters(model, PhaseViscosities = [1e-3, 5e-3]) # 1 and 5 cP
    state0 = setup_state(model, Pressure = p0, Saturations = [0.7, 0.3])
    return (model, state0, parameters, forces, tstep)
end
# Test sensitivity of integrated Buckley-Leverett outlet saturation
G = (model, state, dt, step_no, forces) -> dt*state[:Saturations][2, end]
model, state0, parameters, forces, tstep = setup_bl(nc = 10, nstep = 10)
adj, num = get_sens(model, state0, parameters, tstep, forces, G)

@testset "BL sensitivities" begin
    t = :Transmissibilities
    @testset "$t" begin
        @test isapprox(adj[t], num[t], rtol = 1e-3)
    end
    t = :TwoPointGravityDifference
    @testset "$t" begin
        @test isapprox(adj[t], num[t], rtol = 1e-2)
    end
    t = :FluidVolume
    @testset "$t" begin
        @test isapprox(adj[t], num[t], rtol = 1e-3)
    end
end
