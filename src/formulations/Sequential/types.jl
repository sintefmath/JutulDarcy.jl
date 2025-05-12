struct SequentialSimulator{M, P, T, S} <: Jutul.JutulSimulator
    model::M
    pressure::P
    transport::T
    storage::S
end

function SequentialSimulator(case::JutulCase; kwarg...)
    return SequentialSimulator(case.model; state0 = case.state0, parameters = case.parameters)
end

function SequentialSimulator(model; state0 = setup_state(model), parameters = setup_parameters(model), avg_mobility = false)
    rmodel = reservoir_model(model)
    sys = rmodel.system
    pmodel = convert_to_sequential(model, pressure = true, avg_mobility = avg_mobility)
    tmodel = convert_to_sequential(model, pressure = false)
    function add_total_saturation!(m, state0)
        if !haskey(state0, :TotalSaturation)
            state0[:TotalSaturation] = ones(number_of_cells(m.domain))
        end
    end
    function add_total_saturation!(m::MultiModel, state0)
        add_total_saturation!(m[:Reservoir], state0[:Reservoir])
    end
    add_total_saturation!(model, state0)

    function merge_initial_state(m, state0, parameters)
        return merge(state0, parameters)
    end
    function merge_initial_state(m::MultiModel, state0, parameters)
        init = copy(state0)
        init[:Reservoir] = merge(state0[:Reservoir], parameters[:Reservoir])
        return init
    end
    init = merge_initial_state(model, state0, parameters)
    function subsimulator(m)
        s0, prm = setup_state_and_parameters(m, init)
        return Simulator(m, state0 = s0, parameters = prm)
    end

    PSim = subsimulator(pmodel)
    TSim = subsimulator(tmodel)


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
    S = JutulStorage()
    S[:transfer] = JutulStorage()
    setup_sequential_storage!(S[:transfer], PSim.model, TSim.model)
    S[:recorder] = ProgressRecorder()
    # S[:state] = seq_state
    # S[:state0] = deepcopy(seq_state)
    # S[:primary_variables] = seq_state
    return SequentialSimulator(model, PSim, TSim, S)
end
