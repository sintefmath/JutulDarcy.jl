function JutulDarcy.setup_reservoir_model(afi::AFIInputFile;
        reservoir = JutulDarcy.reservoir_domain(afi),
        system = JutulDarcy.AFISetup.setup_system(afi),
        wells = JutulDarcy.AFISetup.setup_wells(afi, reservoir),
        extra_out = true,
        kwarg...
    )

    model = setup_reservoir_model(reservoir, system; wells = wells, extra_out = false, kwarg...)
    pvars = JutulDarcy.AFISetup.setup_pvt_variables(afi, system, reservoir)
    svars = JutulDarcy.AFISetup.setup_saturation_variables(afi, system, reservoir)
    Jutul.replace_variables!(model; pairs(pvars)..., pairs(svars)..., throw = false)
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
