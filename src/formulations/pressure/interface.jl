function convert_to_sequential(model, variant = :pressure)
    if variant == :pressure
        f = PressureFormulation()
    else
        error("NotImplementedYet")
    end
    # TODO: Figure out context
    pmodel = SimulationModel(
        model.domain,
        model.system,
        data_domain = model.data_domain,
        formulation = f
        )
    for (pkey, pvar) in model.primary_variables
        if pkey != :Pressure
            pmodel.parameters[pkey] = pvar
        end
    end
    for (skey, svar) in model.secondary_variables
        pmodel.secondary_variables[skey] = svar
    end
    for (pkey, pvar) in model.parameters
        pmodel.parameters[pkey] = pvar
    end
    return pmodel
end
