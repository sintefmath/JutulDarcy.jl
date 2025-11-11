function JutulDarcy.setup_reservoir_model(afi::AFIInputFile;
        reservoir = JutulDarcy.reservoir_domain(afi),
        phases = missing,
        system = setup_system(afi, phases = phases),
        wells = setup_wells(afi, reservoir),
        extra_out = false,
        disable_hysteresis = true,
        disable_endscale = false,
        thermal = thermal_type(afi) == :thermal,
        kwarg...
    )

    model = setup_reservoir_model(reservoir, system;
        wells = wells,
        extra_out = false,
        thermal = thermal,
        kwarg...
    )
    pvars = setup_pvt_variables(afi, system, reservoir)
    svars = setup_saturation_variables(afi, system, reservoir,
        disable_hysteresis = disable_hysteresis,
        disable_endscale = disable_endscale,
    )
    allvars = merge(pvars, svars)
    for (k, v) in allvars
        for (model_key, submodel) in pairs(model.models)
            svars0 = Jutul.get_secondary_variables(submodel)
            prm0 = Jutul.get_parameters(submodel)
            if haskey(svars0, k) || haskey(prm0, k)
                Jutul.delete_variable!(submodel, k)
                v::JutulVariables
                submodel.secondary_variables[k] = v
            end
        end
    end
    if haskey(svars, :CapillaryPressure)
        # This should be added to the reservoir model even it is not already a value, check why it does not work.
        model[:Reservoir].secondary_variables[:CapillaryPressure] = svars[:CapillaryPressure]
    end
    JutulDarcy.set_rock_compressibility!(model, afi)
    JutulDarcy.add_relperm_parameters!(model)
    if extra_out
        retval = (model, Jutul.setup_parameters(model))
    else
        retval = model
    end
    return retval
end

function JutulDarcy.setup_reservoir_state(afi::AFIInputFile, model::MultiModel; kwarg...)
    regs = setup_equilibrium_regions(afi, model)
    return JutulDarcy.setup_reservoir_state(model, regs; kwarg...)
end
