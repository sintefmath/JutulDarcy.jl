struct SequentialSimulator{M, P, T, S} <: Jutul.JutulSimulator
    model::M
    pressure::P
    transport::T
    storage::S
end

struct PressureFormulation <: JutulFormulation

end

const PressureModel = SimulationModel{<:Any, <:Any, PressureFormulation, <:Any}

struct PressureReductionFactors <: Jutul.VectorVariables end

Jutul.values_per_entity(model::PressureModel, ::PressureReductionFactors) = number_of_components(model.system)

struct PressureEquation{T} <: JutulEquation where T<:ConservationLaw
    conservation::T
end

struct PressureEquationTPFAStorage{A, HC, B}
    accumulation::A
    accumulation_symbol::Symbol
    half_face_flux_cells::HC
    buf::B
end

function PressureEquationTPFAStorage(model, eq::PressureEquation; ad = true, kwarg...)
    ceq = eq.conservation
    ceq.flow_discretization::TwoPointPotentialFlowHardCoded
    number_of_equations = 1
    D, ctx = model.domain, model.context
    cell_entity = Cells()
    face_entity = Faces()
    nc = count_active_entities(D, cell_entity, for_variables = false)
    nf = count_active_entities(D, face_entity, for_variables = false)
    nhf = number_of_half_faces(ceq.flow_discretization)
    face_partials = degrees_of_freedom_per_entity(model, face_entity)
    @assert face_partials == 0 "Only supported for cell-centered discretization"
    alloc = (n, entity, n_entities_pos) -> CompactAutoDiffCache(number_of_equations, n, model,
                                                                                entity = entity, n_entities_pos = n_entities_pos, 
                                                                                context = ctx; kwarg...)
    # Accumulation terms
    acc = alloc(nc, cell_entity, nc)
    # Source terms - as sparse matrix
    t_acc = eltype(acc.entries)
    # src = sparse(zeros(0), zeros(0), zeros(t_acc, 0), size(acc.entries)...)
    # Half face fluxes - differentiated with respect to pairs of cells
    hf_cells = alloc(nhf, cell_entity, nhf)
    # # Half face fluxes - differentiated with respect to the faces
    # if face_partials > 0
    #     hf_faces = alloc(nf, face_entity, nhf)
    # else
    #     hf_faces = nothing
    # end
    buf = zeros(t_acc, number_of_components(model.system))
    return PressureEquationTPFAStorage(acc, conserved_symbol(ceq), hf_cells, buf)
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
        for (k, m) in pairs(m.models) 
            if m.system isa MultiPhaseSystem
                add_total_saturation!(m, state0[k])
            end
        end
    end
    add_total_saturation!(model, state0)

    function merge_initial_state(m, state0, parameters)
        return merge(state0, parameters)
    end
    function merge_initial_state(m::MultiModel, state0, parameters)
        init = copy(state0)
        for (k, m) in pairs(m.models) 
            if m.system isa MultiPhaseSystem
                init[k] = merge(state0[k], parameters[k])
            end
        end
        return init
    end
    init = merge_initial_state(model, state0, parameters)
    function subsimulator(m)
        s0, prm = setup_state_and_parameters(m, init)
        return Simulator(m, state0 = s0, parameters = prm)
    end

    reservoir_model(pmodel)::PressureModel
    reservoir_model(tmodel)::TransportModel
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
