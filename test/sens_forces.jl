using Jutul, JutulDarcy, Test
# import SparsityTracing as ST
function setup_bl_twoforces(;nc = 100, time = 1.0, nstep = 100)
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
    src  = SourceTerm(1, irate, fractional_flow = [1.0, 0.0])
    bc = FlowBoundaryCondition(nc, p0/2)
    forces = setup_forces(model, sources = src, bc = bc)

    src2 = SourceTerm(1, 0.99*irate, fractional_flow = [1.0, 0.0])
    force2 = setup_forces(model, sources = src2, bc = bc)
    forces = repeat([forces], length(tstep))
    for i in 1:(nstep÷2)
        forces[i] = force2
    end

    parameters = setup_parameters(model, PhaseViscosities = [1e-3, 5e-3]) # 1 and 5 cP
    state0 = setup_state(model, Pressure = p0, Saturations = [0.7, 0.3])
    return (model, state0, parameters, forces, tstep)
end

function numerical_diff_bl(model, state0, parameters, forces, tstep, G)
    dx = Float64[]
    s0, = simulate(state0, model, tstep, forces = forces, parameters = parameters, info_level = -1)
    obj0 = Jutul.evaluate_objective(G, model, s0, tstep, forces)
    ϵ = 1e-6
    unique_forces, to_step = Jutul.unique_forces_and_mapping(forces, tstep)
    for fno in eachindex(unique_forces)
        x, cfg = Jutul.vectorize_forces(unique_forces[fno], model)
        for i in eachindex(x)
            x_delta = copy(x)
            x_delta[i] += ϵ
            new_force = Jutul.devectorize_forces(unique_forces[fno], model, x_delta, cfg)
            new_forces = deepcopy(forces)
            for j in to_step[fno]
                new_forces[j] = new_force
            end
            s, r = simulate(state0, model, tstep, forces = new_forces, parameters = parameters, info_level = -1)
            obj = Jutul.evaluate_objective(G, model, s, tstep, new_forces)
            push!(dx, (obj - obj0)/ϵ)
        end
    end
    return dx
end


model, state0, parameters, forces, tstep = setup_bl_twoforces(nc = 10, nstep = 10)
G = (model, state, dt, step_no, forces) -> dt*(sum(state[:Saturations][1, :] .- 0.5))^2

dx = numerical_diff_bl(model, state0, parameters, forces, tstep, G)

states, reports = simulate(state0, model, tstep, forces = forces, parameters = parameters)

dforces, grad_adj = Jutul.solve_adjoint_forces(model, states, reports, G, forces,
                state0 = state0, parameters = parameters)


# case = JutulCase(model, tstep, forces, state0 = state0, parameters = parameters)
## opt_config = Jutul.forces_optimization_config(model, forces, tstep, :all, abs_min = 0.0)
# x0, xmin, xmax, f, g!, out = Jutul.setup_force_optimization(case, G, opt_config)
