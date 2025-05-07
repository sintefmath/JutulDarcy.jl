
function convert_to_sequential(model; avg_mobility = false, pressure = true)
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
    for (skey, svar) in model.secondary_variables
        seqmodel.secondary_variables[skey] = svar
    end
    if pressure
        for (pkey, pvar) in model.primary_variables
            if pkey != :Pressure
                seqmodel.parameters[pkey] = pvar
            end
        end
        if avg_mobility
            mob = :PhaseMobilities
            if haskey(seqmodel.secondary_variables, mob)
                seqmodel.parameters[mob] = PhaseMobilities()
                delete!(seqmodel.secondary_variables, mob)
            end
        end
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
    for (pkey, pvar) in seqmodel.parameters
        if haskey(seqmodel.secondary_variables, pkey)
            delete!(seqmodel.parameters, pkey)
        end
    end
    return seqmodel
end


function convert_to_sequential(model::MultiModel; pressure = true, kwarg...)
    ct = Vector{Jutul.CrossTermPair}()
    for ctp in model.cross_terms
        ctp = deepcopy(ctp)
        (; target, source, target_equation, source_equation, cross_term) = ctp
        target_is_cons = target_equation == :mass_conservation
        source_is_cons = source_equation == :mass_conservation
        if pressure && (target_is_cons || source_is_cons)
            if source_is_cons
                source_equation = :pressure
            end
            if target_is_cons
                target_equation = :pressure
            end
            if target_is_cons && source_is_cons && cross_term isa ReservoirFromWellFlowCT
                cross_term = PressureReservoirFromWellFlowCT(cross_term)
            end
            ctp = Jutul.CrossTermPair(target, source, target_equation, source_equation, cross_term)
        end
        push!(ct, ctp)
    end
    smodel = convert_to_sequential(model[:Reservoir]; pressure = pressure, kwarg...)
    models = Dict{Symbol, Any}()
    for (k, v) in pairs(model.models)
        if k == :Reservoir
            models[k] = smodel
        else
            if v.system isa MultiPhaseSystem && pressure
                v = convert_to_sequential(v; pressure = true, kwarg...)
            end
            models[k] = deepcopy(v)
        end
    end

    if isnothing(model.groups)
        g = nothing
    else
        g = copy(model.groups)
    end
    seqmodel = MultiModel(
        models,
        cross_terms = ct,
        groups = g,
        context = model.context,
        reduction = model.reduction
        )
    return seqmodel
end
