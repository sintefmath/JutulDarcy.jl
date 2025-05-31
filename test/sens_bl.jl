using Jutul, JutulDarcy, Test

function get_sens(model, state0, parameters, tstep, forces, G)
    sim = Simulator(model, state0 = state0, parameters = parameters)
    states, reports = simulate(sim, tstep, forces = forces, extra_timing = false, info_level = -1)

    grad_adj = solve_adjoint_sensitivities(model, states, reports, G,
                    forces = forces, state0 = state0, parameters = parameters,
                    extra_timing = false, raw_output = false)

    grad_numeric = Dict{Symbol, Any}()

    is_multi = isa(model, MultiModel)
    if is_multi
        grad_adj = grad_adj[:Reservoir]
    end
    for k in keys(grad_adj)
        if is_multi
            target = (:Reservoir, k)
        else
            target = k
        end
        grad_numeric[k] = Jutul.solve_numerical_sensitivities(model, states, reports, G, target,
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
    kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0], [0.2, 0.2])
    replace_variables!(model, RelativePermeabilities = kr)
    tot_time = sum(timesteps)
    pv = pore_volume(G)
    irate = 500*sum(pv)/tot_time
    src  = [SourceTerm(1, irate, fractional_flow = [0.8, 0.2]), 
            SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
    forces = setup_forces(model, sources = src)

    parameters = setup_parameters(model, PhaseViscosities = [1e-3, 5e-3]) # 1 and 5 cP
    state0 = setup_state(model, Pressure = p0, Saturations = [0.7, 0.3])
    return (model, state0, parameters, forces, tstep)
end
# Test sensitivity of integrated Buckley-Leverett outlet saturation
@testset "BL sensitivities" begin
    for model_type in [:single, :multi, :multi_spec, :multi_dummy]
        model, state0, parameters, forces, tstep = setup_bl(nc = 10, nstep = 10)
        if model_type == :multi_dummy
            model = MultiModel(Dict(:Dummy => model, :Reservoir => model))
            forces = Dict(:Dummy => forces, :Reservoir => forces)
            parameters = Dict(:Dummy => parameters, :Reservoir => parameters)
            state0 = Dict(:Dummy => state0, :Reservoir => state0)
            G = (m, state, dt, step_info, forces) -> dt*state[:Reservoir][:Saturations][2, end]
        elseif model_type == :multi || model_type == :multi_spec
            model = MultiModel(Dict(:Reservoir => model), specialize = model_type == :multi_spec)
            forces = Dict(:Reservoir => forces)
            parameters = Dict(:Reservoir => parameters)
            state0 = Dict(:Reservoir => state0)
            G = (m, state, dt, step_info, forces) -> dt*state[:Reservoir][:Saturations][2, end]
        else
            @assert model_type == :single
            G = (model, state, dt, step_info, forces) -> dt*state[:Saturations][2, end]
        end
        adj, num = get_sens(model, state0, parameters, tstep, forces, G)
        @testset "$model_type" begin
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
    end
end
