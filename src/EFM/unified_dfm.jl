function setup_unified_dfm_model(matrix, fractures, system;
        wells = [], extend_wells = true, subtract_fracture_volume = true,
        extra_out = false, kwarg...)
    # Initialize the same coefficients and well completions as coupled DFM.
    coupled = setup_fractured_reservoir_model(copy_fracture_data(matrix),
        copy_fracture_data(fractures), system;
        wells = [copy_fracture_data(well) for well in wells],
        extend_wells = extend_wells, subtract_fracture_volume = subtract_fracture_volume,
        kwarg...)
    parameters = setup_parameters(coupled)
    domain = unified_dfm_domain(coupled, parameters)
    nm = number_of_cells(matrix)
    wells = [unified_fracture_well(well, parameters[name], nm)
        for (name, well) in pairs(get_model_wells(coupled; data_domain = true))]
    model = setup_reservoir_model(domain, system; wells = wells, kwarg...)
    if JutulDarcy.model_is_thermal(reservoir_model(model))
        set_parameters!(reservoir_model(model), FluidThermalConductivities = EFMFluidThermalConductivities())
    end
    out = model
    if extra_out
        out = (model, setup_parameters(model))
    end
    return out
end

function unified_dfm_domain(model, parameters)
    matrix = reservoir_domain(model)
    fractures = model[:Fractures].data_domain
    fracture_mesh = physical_representation(fractures)
    fmc = FractureMatrixConnection()
    mc, fc = fractures[:matrix_cells, fmc], fractures[:connection_cells, fmc]
    mf = fractures[:matrix_faces, fmc]
    areas = tpfv_geometry(fracture_mesh).volumes[fc]
    domain = unified_fracture_domain(matrix, fractures, areas, matrix[:normals][:, mf])
    pm, pf = parameters[:Reservoir], parameters[:Fractures]
    domain[:transmissibilities, Faces()] = vcat(pm[:Transmissibilities],
        pf[:Transmissibilities], pf[:FractureMatrixTransmissibility])
    if JutulDarcy.model_is_thermal(reservoir_model(model))
        # DFM cross terms use total heat conductivity independent of saturation.
        # Store it in the rock term, with zero fluid conductivity on these faces.
        domain[:rock_thermal_conductivities, Faces()] = vcat(pm[:RockThermalConductivities],
            pf[:RockThermalConductivities], pf[:FractureMatrixThermalConductivity])
        nph = size(pm[:FluidThermalConductivities], 1)
        domain[:fluid_thermal_conductivities, Faces()] = hcat(pm[:FluidThermalConductivities],
            pf[:FluidThermalConductivities], zeros(eltype(domain[:volumes]), nph, length(mc)))
    end
    return domain
end
