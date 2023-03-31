
function convert_to_sequential(model, variant = :pressure)
    is_pressure = variant == :pressure
    if is_pressure
        f = PressureFormulation()
    else
        @assert variant == :transport
        f = TransportFormulation()
    end
    # TODO: Figure out context
    pmodel = SimulationModel(
        model.domain,
        model.system,
        data_domain = model.data_domain,
        formulation = f
        )
    if is_pressure
        for (pkey, pvar) in model.primary_variables
            if pkey != :Pressure
                pmodel.parameters[pkey] = pvar
            end
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
