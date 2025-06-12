function setup_reservoir_dict_optimization(dict, F = missing)
    return DictParameters(dict, F)
end

function setup_reservoir_dict_optimization(case::JutulCase;
        use_trans = false,
        use_pore_volume = false,
        strict = false,
        kwarg...
    )
    is_multimodel = case.model isa MultiModel
    if is_multimodel
        rparameters = case.parameters[:Reservoir]
        state0_r = case.state0[:Reservoir]
    else
        rparameters = case.parameters
        state0_r = case.state0
    end
    rmodel = reservoir_model(case.model)
    rdomain = reservoir_domain(rmodel)

    skip_list = [
        :satnum,
        :eqlnum,
        :pvtnum,
        :cell_centroids,
        :volumes,
        :neighbors,
        :areas,
        :normals,
        :face_centroids,
        :half_face_cells,
        :half_face_faces,
        :boundary_areas,
        :boundary_centroids,
        :boundary_normals,
        :boundary_neighbors,
    ]
    is_thermal = model_is_thermal(case.model)
    if !is_thermal
        push!(skip_list, :rock_heat_capacity)
        push!(skip_list, :component_heat_capacity)
        push!(skip_list, :rock_heat_capacity)
        push!(skip_list, :fluid_thermal_conductivity)
        push!(skip_list, :rock_thermal_conductivity)
    end
    skip_list_parameters = Symbol[:Transmissibilities]
    DT = OrderedDict{Symbol, Any}
    opt_dict = DT()
    prm_dict = DT()
    dd_dict = DT()
    state0_dict = DT()
    w_dict = DT()
    opt_dict[:model] = dd_dict
    opt_dict[:parameters] = prm_dict
    opt_dict[:state0] = state0_dict
    opt_dict[:wells] = w_dict
    if use_trans
        opt_dict[:Transmissibilities] = copy(rparameters[:Transmissibilities])
        push!(skip_list, :permeability)
    end
    if use_pore_volume
        opt_dict[:pore_volume] = pore_volume(rdomain)
        push!(skip_list, :porosity)
        push!(skip_list, :net_to_gross)
    end
    for k in setdiff(keys(rdomain), skip_list)
        v = rdomain[k]
        if eltype(v)<:AbstractFloat
            dd_dict[k] = v
        end
    end
    if is_multimodel
        for (k, submodel) in pairs(case.model.models)
            if model_or_domain_is_well(submodel)
                subdict = DT()
                subprm = case.parameters[k]
                subdict[:WellIndices] = subprm[:WellIndices]
                if haskey(subprm, :WellIndicesThermal) && is_thermal
                    subdict[:WellIndicesThermal] = subprm[:WellIndicesThermal]
                end
                w_dict[k] = subdict
            end
        end
    end
    # Initial conditions...
    # TODO: Write F
    F = missing
    return DictParameters(opt_dict, F, strict = strict)
end

function optimize_reservoir(dopt, objective, setup_fn = dopt.setup_function;
        simulator_arg = (output_substates = false, ),
        kwarg...
    )
    if ismissing(setup_fn)
        error("Setup function was not found in DictParameters struct or as last positional argument.")
    end
    case0 = setup_fn(dopt.parameters, missing)
    sim, cfg = setup_reservoir_simulator(case0; info_level = -1, simulator_arg...)
    return Jutul.optimize(dopt, objective, setup_fn; simulator = sim, config = cfg, kwarg...)
end

function parameters_gradient_reservoir(dopt, objective, setup_fn = dopt.setup_function;
        simulator_arg = (output_substates = false, ),
        kwarg...
    )
    if ismissing(setup_fn)
        error("Setup function was not found in DictParameters struct or as last positional argument.")
    end
    case0 = setup_fn(dopt.parameters, missing)
    sim, cfg = setup_reservoir_simulator(case0; info_level = -1, simulator_arg...)
    return Jutul.parameters_gradient(dopt, objective, setup_fn; simulator = sim, config = cfg, kwarg...)
end
