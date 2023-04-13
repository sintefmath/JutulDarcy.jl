
function convert_to_sequential(model; pressure = true)
    if pressure
        f = PressureFormulation()
    else
        f = TransportFormulation()
    end
    transport = !pressure
    # TODO: Figure out context
    seqmodel = SimulationModel(
        model.domain,
        model.system,
        data_domain = model.data_domain,
        formulation = f
        )
    if pressure
        for (pkey, pvar) in model.primary_variables
            if pkey != :Pressure
                seqmodel.parameters[pkey] = pvar
            end
        end
    end
    for (skey, svar) in model.secondary_variables
        seqmodel.secondary_variables[skey] = svar
    end
    if transport
        vars = seqmodel.secondary_variables
        if haskey(vars, :ShrinkageFactors)
            k = :ShrinkageFactors
        else
            k = :PhaseMassDensities
        end
        @assert haskey(vars, k)
        vars[:UncorrectedVariable] = vars[k]
        vars[k] = TotalSaturationCorrectedVariable()
        push!(seqmodel.output_variables, :Pressure)
    end
    for (pkey, pvar) in model.parameters
        seqmodel.parameters[pkey] = pvar
    end
    return seqmodel
end


function convert_to_sequential(model::MultiModel; kwarg...)
    ct = deepcopy(model.cross_terms)
    smodel = convert_to_sequential(model[:Reservoir]; kwarg...)
    models = Dict{Symbol, Any}()
    for (k, v) in pairs(model.models)
        if k == :Reservoir
            models[k] = smodel
        else
            models[k] = deepcopy(v)
        end
    end

    seqmodel = MultiModel(
        models,
        cross_terms = ct,
        groups = copy(model.groups),
        context = model.context,
        reduction = model.reduction
        )
    return seqmodel
end
