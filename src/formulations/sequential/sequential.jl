struct SequentialSimulator{M, P, T, S} <: Jutul.JutulSimulator
    model::M
    pressure::P
    transport::T
    storage::S
end

function SequentialSimulator(model; state0 = setup_state(model), parameters = setup_parameters(model))
    pmodel = convert_to_sequential(model, pressure = true)
    tmodel = convert_to_sequential(model, pressure = false)
    init = merge(state0, parameters)
    if !haskey(state0, :TotalSaturation)
        init[:TotalSaturation] = ones(number_of_cells(model.domain))
    end
    function subsimulator(m)
        s0, prm = setup_state_and_parameters(m; pairs(init)...)
        return Simulator(m, state0 = s0, parameters = prm)    
    end

    PSim = subsimulator(pmodel)
    TSim = subsimulator(tmodel)
    S = JutulStorage()
    S[:recorder] = ProgressRecorder()
    return SequentialSimulator(model, PSim, TSim, S)
end

function Jutul.select_linear_solver(sim::SequentialSimulator; kwarg...)
    return nothing
end

# function Jutul.reset_variables!(sim::SequentialSimulator, vars; kwarg...)
#     Jutul.reset_variables!(sim.pressure, vars; kwarg...)
#     Jutul.reset_variables!(sim.transport, vars; kwarg...)
#     return sim
# end

function Jutul.initial_setup!(sim::SequentialSimulator, config, timesteps; kwarg...)
    Jutul.initial_setup!(sim.pressure, config, timesteps; kwarg...)
end

function Jutul.initialize_before_first_timestep!(sim::SequentialSimulator, dt; kwarg...)
    Jutul.initialize_before_first_timestep!(sim.pressure, dt; kwarg...)
end

function Jutul.perform_step!(
        simulator::SequentialSimulator,
        dt,
        forces,
        config;
        iteration = NaN,
        relaxation = 1.0,
        update_secondary = true,
        solve = true
    )
    # Solve pressure
    max_iter = config[:max_nonlinear_iterations]
    done_p, report_p = Jutul.solve_ministep(simulator.pressure, dt, forces, max_iter, config)
    if done_p
        # Copy over values for pressure and fluxes into parameters for second simulator

        # Then transport
        done_t, report_t = Jutul.solve_ministep(simulator.transport, dt, forces, max_iter, config)
    else
        error("Pressure failure not implemented")
    end
    # Return convergence criterion for outer loop
end


include("interface.jl")
