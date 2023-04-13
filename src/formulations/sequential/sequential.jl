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

    function seq_output_keys(m)
        return copy(m.output_variables)
    end
    # Transfer rule: All parameters present in other that are not parameters both places.
    function transfer_keys(target, source)
        prm_t = keys(Jutul.get_parameters(target))
        prm_s = keys(Jutul.get_parameters(source))

        all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        return intersect(all_s, setdiff(prm_t, prm_s))
    end
    transfer_keys_transport = transfer_keys(TSim.model, PSim.model)
    transfer_keys_pressure = transfer_keys(PSim.model, TSim.model)
    # Init rule: Transfer over everything that is present from transport state.
    function init_keys(target, source)
        all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        return intersect(all_t, all_s)
    end
    # init_keys_transport = init_keys(TSim.model)
    init_keys_pressure = init_keys(PSim.model, TSim.model)

    transfer_keys = JutulStorage()
    transfer_keys[:pressure] = transfer_keys_pressure
    transfer_keys[:transport] = transfer_keys_transport

    init_keys = JutulStorage()
    init_keys[:pressure] = init_keys_pressure
    # init_keys[:transport] = init_keys_transport

    S = JutulStorage()
    S[:init_keys] = init_keys
    S[:transfer_keys] = transfer_keys

    # @info "Keys" transfer_keys_pressure transfer_keys_transport init_keys_transport init_keys_pressure


    # all_keys = union(seq_output_keys(PSim.model), seq_output_keys(TSim.model))
    # seq_state = JutulStorage()
    # for s in [PSim, TSim]
    #     s0 = s.storage.state0
    #     for k in seq_output_keys(s.model)
    #         if haskey(s0, k)
    #             seq_state[k] = similar(s0[k])
    #         end
    #     end
    # end
    # @info "!" seq_state
    # error()

    S[:recorder] = ProgressRecorder()
    # S[:state] = seq_state
    # S[:state0] = deepcopy(seq_state)
    # S[:primary_variables] = seq_state
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
    Jutul.initial_setup!(sim.transport, config, timesteps; kwarg...)
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
    psim = simulator.pressure
    pstate = psim.storage.state

    tsim = simulator.transport
    tstate = tsim.storage.state
    
    transfer_keys = simulator.storage.transfer_keys
    # Solve pressure
    max_iter = config[:max_nonlinear_iterations]
    config[:always_update_secondary] = true
    # Copy over variables to parameters for both solves
    if iteration > 1
        # We need to transfer from pressure
        for k in transfer_keys[:pressure]
            update_values!(pstate[k], tstate[k])
        end
    end
    done_p, report_p = Jutul.solve_ministep(psim, dt, forces, max_iter, config)
    if done_p
        # Copy over values for pressure and fluxes into parameters for second simulator
        model_p = psim.model
        state_p = psim.storage.state
        vT = tsim.storage.state.TotalVolumetricFlux
        store_total_fluxes!(vT, model_p, as_value(state_p))
        for k in transfer_keys[:transport]
            update_values!(tstate[k], pstate[k])
        end
        # Then transport
        done_t, report_t = Jutul.solve_ministep(tsim, dt, forces, max_iter, config)
    else
        error("Pressure failure not implemented")
    end
    # Return convergence criterion for outer loop if SFI
    return (0.0, done_t, report_p)
end

function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    report = Dict{Symbol, Any}()
    report[:pressure] = Jutul.update_after_step!(sim.pressure, dt, forces; kwarg...)
    report[:transport] = Jutul.update_after_step!(sim.transport, dt, forces; kwarg...)
    return report
end

function Jutul.get_output_state(sim::SequentialSimulator)
    Jutul.get_output_state(sim.transport)
end

function Jutul.final_simulation_message(simulator::SequentialSimulator, arg...)
    Jutul.jutul_message("Sequential", "All done.")
    # Jutul.final_simulation_message(simulator.pressure, arg...)
    # Jutul.jutul_message("Transport")
    #Jutul.final_simulation_message(simulator.transport, arg...)
end

function Jutul.update_before_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.pressure, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.transport, dt, forces; kwarg...)
end

function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    # First, sync up state and state0 for transport
    Jutul.update_after_step!(sim.transport, dt, forces; kwarg...)
    # NOTE: This part must be done carefully so that the pressure contains the
    # final solution from the last transport solve.
    pstate = sim.pressure.storage.state
    tstate = sim.transport.storage.state
    for k in sim.storage.init_keys.pressure
        update_values!(pstate[k], tstate[k])
    end
    Jutul.update_after_step!(sim.pressure, dt, forces; kwarg...)
end


function Jutul.reset_state_to_previous_state!(sim::SequentialSimulator)
    Jutul.reset_state_to_previous_state!(sim.pressure)
    Jutul.reset_state_to_previous_state!(sim.transport)
end


function Jutul.reset_variables!(sim::SequentialSimulator, state; kwarg...)
    @info "Hey..."
    Jutul.reset_variables!(sim.pressure, state; kwarg...)
    Jutul.reset_variables!(sim.transport, state; kwarg...)
end

function Jutul.reset_previous_state!(sim::SequentialSimulator, state)
    Jutul.reset_previous_state!(sim.pressure, state)
    Jutul.reset_previous_state!(sim.transport, state)
end
# reset_previous_state!

include("interface.jl")
