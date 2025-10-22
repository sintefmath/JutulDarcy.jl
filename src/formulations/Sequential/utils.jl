

function setup_sequential_storage!(S, p_model::MultiModel, t_model::MultiModel)
    for k in Jutul.submodels_symbols(p_model)
        p_submodel = p_model[k]
        t_submodel = t_model[k]
        is_reservoir = k == :Reservoir
        S[k] = JutulStorage()
        setup_sequential_storage!(S[k], p_submodel, t_submodel,
            is_reservoir = is_reservoir,
            is_well = model_or_domain_is_well(p_submodel)
        )
    end
end

function setup_sequential_storage!(S, p_model::SimulationModel, t_model::SimulationModel;
        is_reservoir = true,
        is_well = !is_reservoir
    )
    function seq_output_keys(m)
        return copy(m.output_variables)
    end
    # Transfer rule: All parameters present in other that are not parameters
    # both places and all primary variables that are primary in both models.
    function transfer_keys(target, source)
        prm_t = keys(Jutul.get_parameters(target))
        prm_s = keys(Jutul.get_parameters(source))

        # all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        parameter_overlap = intersect(all_s, setdiff(prm_t, prm_s))

        pvar_t = keys(Jutul.get_variables_by_type(target, :primary))
        pvar_s = keys(Jutul.get_variables_by_type(source, :primary))

        primary_overlap = intersect(pvar_t, pvar_s)
        return [parameter_overlap..., primary_overlap...]
    end
    transfer_keys_transport = transfer_keys(t_model, p_model)
    transfer_keys_pressure = transfer_keys(p_model, t_model)
    # Init rule: Transfer over everything that is present from transport state.
    function init_keys(target, source)
        all_t = keys(Jutul.get_variables_by_type(target, :all))
        all_s = keys(Jutul.get_variables_by_type(source, :all))

        # return intersect(all_t, all_s)

        secondary_s = keys(Jutul.get_variables_by_type(source, :secondary))
        tkeys = setdiff(intersect(all_t, all_s), secondary_s)
        if is_well || is_reservoir
            push!(tkeys, :TotalMasses)
        end
        return tkeys
    end
    init_keys_pressure = init_keys(p_model, t_model)

    transfer_keys = JutulStorage()
    transfer_keys[:pressure] = transfer_keys_pressure
    transfer_keys[:transport] = transfer_keys_transport

    init_keys = JutulStorage()
    init_keys[:pressure] = init_keys_pressure
    # init_keys[:transport] = init_keys_transport

    S[:init_keys] = init_keys
    S[:transfer_keys] = transfer_keys

    if is_reservoir
        nph = number_of_phases(t_model.system)
        nc = number_of_cells(t_model.domain)
        λ = zeros(nph, nc)
        @. λ = NaN
        S[:mobility] = λ
        S[:mobility_prev] = similar(λ)
    end
    return S
end

function Jutul.simulator_config(
        sim::SequentialSimulator;
        transport_substeps = 1,
        saturation_tol = 1e-2,
        mobility_tol = 1e-2,
        sfi = false,
        linear_solver = missing,
        kwarg...
    )
    cfg = Jutul.JutulConfig("Simulator config")
    Jutul.add_option!(cfg, :transport_substeps, transport_substeps, types = Int)
    Jutul.add_option!(cfg, :saturation_tol, saturation_tol, types = Float64)
    Jutul.add_option!(cfg, :mobility_tol, mobility_tol, types = Float64)
    Jutul.add_option!(cfg, :sfi, sfi, types = Bool)

    Jutul.simulator_config!(cfg, sim;
        min_nonlinear_iterations = 0,
        kwarg...
    )

    for k in [:pressure, :transport]
        cfg_k = Jutul.simulator_config(getfield(sim, k); always_update_secondary = true, kwarg...)
        Jutul.add_option!(cfg, k, cfg_k)
    end
    return cfg
end

function Jutul.simulator_executor(sim::SequentialSimulator)
    return Jutul.simulator_executor(sim.pressure)
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

function get_reservoir_state(sim, current = true)
    if current
        model_state = sim.storage.state
    else
        model_state = sim.storage.state0
    end
    if sim.model isa MultiModel
        model_state = model_state.Reservoir
    end
    return model_state
end


function sequential_sync_values!(sim::SequentialSimulator; to_key::Symbol = :pressure, kind = :transfer_keys)
    @assert to_key in (:pressure, :transport)
    pstate = sim.pressure.storage.state
    tstate = sim.transport.storage.state
    if to_key == :pressure
        dest = pstate
        src = tstate
        model = sim.transport.model
    else
        dest = tstate
        src = pstate
        model = sim.pressure.model
    end
    sequential_sync_values!(dest, src, model, sim.storage.transfer, to_key, kind)
end

function sequential_sync_values!(dest, src, model::SimulationModel, storage_transfer, dest_key, kind)
    for k in storage_transfer[kind][dest_key]
        # @info "Syncing $k from $dest_key" src[k]
        update_values!(dest[k], src[k])
    end
end

function sequential_sync_values!(dest, src, model::MultiModel, storage_transfer, k, kind)
    for modname in keys(storage_transfer)
        sequential_sync_values!(dest[modname], src[modname], model[modname], storage_transfer[modname], k, kind)
    end
end

# function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
#     report = Dict{Symbol, Any}()
#     report[:pressure] = Jutul.update_after_step!(sim.pressure, dt, forces; kwarg...)
#     report[:transport] = Jutul.update_after_step!(sim.transport, dt, forces; kwarg...)
#     return report
# end

function Jutul.get_output_state(sim::SequentialSimulator)
    Jutul.get_output_state(sim.transport)
end

function Jutul.final_simulation_message(simulator::SequentialSimulator, p, rec, t_elapsed, reports, timesteps, config, arg...)
    info_level = config[:info_level]
    print_end_report = config[:end_report]
    verbose = info_level >= 0
    do_print = verbose || print_end_report
    if do_print
        Jutul.jutul_message("Sequential", "Total timing, per SFI iteration")
        Jutul.final_simulation_message(simulator.pressure, p, rec, t_elapsed, reports, timesteps, config, arg...)
        stats = Dict()
        stats[:pressure] = Dict(:iterations => 0, :linear_iterations => 0)
        stats[:transport] = Dict(:iterations => 0, :linear_iterations => 0)

        for rep in reports
            for ministep in rep[:ministeps]
                for step in ministep[:steps]
                    for k in [:pressure, :transport]
                        if !haskey(step, k)
                            continue
                        end
                        for sstep in step[k][:steps]
                            if haskey(sstep, :update)
                                stats[k][:iterations] += 1
                                stats[k][:linear_iterations] += sstep[:linear_iterations]
                            end
                        end
                    end
                end
            end
        end
        for k in [:pressure, :transport]
            s = stats[k]
            i = s[:iterations]
            li = s[:linear_iterations]
            Jutul.jutul_message(titlecase("$k"), "$i iterations, $li linear iterations")
        end
    end
end

function Jutul.update_before_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.pressure, dt, forces; kwarg...)
    Jutul.update_before_step!(sim.transport, dt, forces; kwarg...)
    sequential_set_state0_values!(sim.storage.state0, sim.transport.model, sim.transport.storage.state0)
end

function sequential_set_state0_values!(state_sim, m::SimulationModel, state)
    pvar = Jutul.get_primary_variables(m)
    for k in keys(pvar)
        state_sim[k] .= state[k]
    end
    svar = Jutul.get_secondary_variables(m)
    for k in keys(svar)
        state_sim[k] .= state[k]
    end
    return state_sim
end

function sequential_set_state0_values!(state_sim, m::MultiModel, state)
    for k in Jutul.submodels_symbols(m)
        sequential_set_state0_values!(state_sim[k], m[k], state[k])
    end
    return state_sim
end

function Jutul.update_after_step!(sim::SequentialSimulator, dt, forces; kwarg...)
    # NOTE: This part must be done carefully so that the pressure contains the
    # final solution from the last transport solve.
    sequential_sync_values!(sim, to_key = :pressure, kind = :init_keys)
    return sequential_update_report(sim.transport.model, sim.transport.storage.state, sim.storage.state0)
end

function sequential_update_report(m::SimulationModel, state, state0)
    report = Jutul.JUTUL_OUTPUT_TYPE()
    pvar = Jutul.get_primary_variables(m)
    for k in keys(pvar)
        report[k] = Jutul.variable_change_report(state[k], state0[k], pvar[k])
    end
    svar = Jutul.get_secondary_variables(m)
    for k in keys(svar)
        report[k] = Jutul.variable_change_report(state[k], state0[k], svar[k])
    end
    return report
end

function sequential_update_report(m::MultiModel, state, state0)
    report = Jutul.JUTUL_OUTPUT_TYPE()
    for k in Jutul.submodels_symbols(m)
        report[k] = sequential_update_report(m[k], state[k], state0[k])
    end
    return report
end


function Jutul.reset_state_to_previous_state!(sim::SequentialSimulator)
    Jutul.reset_state_to_previous_state!(sim.pressure)
    Jutul.reset_state_to_previous_state!(sim.transport)
end


function Jutul.reset_variables!(sim::SequentialSimulator, state; kwarg...)
    Jutul.reset_variables!(sim.pressure, state; kwarg...)
    Jutul.reset_variables!(sim.transport, state; kwarg...)
end

function Jutul.reset_previous_state!(sim::SequentialSimulator, state)
    Jutul.reset_previous_state!(sim.pressure, state)
    Jutul.reset_previous_state!(sim.transport, state)
end
# reset_previous_state!

include("interface.jl")
