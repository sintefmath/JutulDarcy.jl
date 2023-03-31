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


function Jutul.initial_setup!(sim::SequentialSimulator, config, timesteps; kwarg...)
    Jutul.initial_setup!(sim.pressure, config, timesteps; kwarg...)
end

function Jutul.initialize_before_first_timestep!(sim::SequentialSimulator, dt; kwarg...)
    Jutul.initialize_before_first_timestep!(sim.pressure, dt; kwarg...)
end

function Jutul.perform_step!(simulator::SequentialSimulator, dt, forces, config; kwarg...)
    # Solve pressure
    e_p, converged_p, report_p = Jutul.perform_step!(simulator.pressure, dt, forces, config; kwarg...)
    # Copy over values for pressure and fluxes into parameters for second simulator
    # Then transport
    # Return convergence criterion for outer loop
end


include("interface.jl")
