function setup_reservoir_dict_optimization(dict, F = missing)
    return DictParameters(dict, F)
end

function setup_reservoir_dict_optimization(case::JutulCase;
        use_trans = false,
        use_pore_volume = false,
        strict = false,
        do_copy = true,
        parameters = Symbol[],
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
    opt_dict[:model] = dd_dict
    opt_dict[:parameters] = prm_dict
    if use_trans
        prm_dict[:Transmissibilities] = copy(rparameters[:Transmissibilities])
        push!(skip_list, :permeability)
    end
    if use_pore_volume
        if haskey(rparameters, :StaticFluidVolume)
            k = :StaticFluidVolume
        else
            k = :FluidVolume
        end
        dd_dict[:pore_volume] = copy(rparameters[k])
        push!(skip_list, :porosity)
        push!(skip_list, :net_to_gross)
        push!(skip_list_parameters, k)
    end
    for k in setdiff(keys(rdomain), skip_list)
        v = rdomain[k]
        if eltype(v)<:AbstractFloat
            dd_dict[k] = copy(v)
        end
    end
    # Regular parameters - requested
    for k in setdiff(parameters, skip_list_parameters)
        if haskey(rparameters, k)
            v = copy(rparameters[k])
            prm_dict[k] = v
        else
            @warn "Parameter $k not found in reservoir model parameters."
        end
    end
    # Wells
    w_dict = DT()
    opt_dict[:wells] = w_dict
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
    # state0
    state0_dict = DT()
    for (k, var) in pairs(rmodel.primary_variables)
        if k == :BlackOilUnknown
            bo_unknown = state0_r[k]
            swat = get(state0_r, :ImmiscibleSaturation, missing)
            sw_fun(x::Missing, i) = 0.0
            sw_fun(x, i) = x[i]
            p = state0_r[:Pressure]
            nc = length(p)
            ix = eachindex(p)
            if has_disgas(rmodel.system)
                rsdef = rmodel.secondary_variables[:Rs]
                rs = similar(p)
                JutulDarcy.update_rs!(rs, rsdef, rmodel, p, bo_unknown, ix)
                state0_dict[:Rs] = rs
            end
            if has_vapoil(rmodel.system)
                rvdef = rmodel.secondary_variables[:Rv]
                rv = similar(p)
                JutulDarcy.update_rv!(rv, rvdef, rmodel, p, bo_unknown, ix)
                state0_dict[:Rv] = rv
            end
            satdef = rmodel.secondary_variables[:Saturations]
            nph = number_of_phases(rmodel.system)
            s = zeros(nph, nc)
            if nph == 3
                JutulDarcy.update_saturations!(s, satdef, rmodel, bo_unknown, state0_r[:ImmiscibleSaturation], ix)
                _, l, v = phase_indices(rmodel.system)
            else
                JutulDarcy.update_saturations!(s, satdef, rmodel, bo_unknown, ix)
                l, v = phase_indices(rmodel.system)
            end
            state0_dict[:LiquidSaturation] = s[l, :]
            state0_dict[:VaporSaturation] = s[v, :]
        else
            state0_dict[k] = copy(state0_r[k])
        end
    end
    opt_dict[:state0] = state0_dict
    F(D, step_info = missing) = optimization_resetup_reservoir_case(D, case, step_info, do_copy = do_copy)
    return DictParameters(opt_dict, F, strict = strict)
end

function optimization_resetup_reservoir_case(opt_dict::AbstractDict, case::JutulCase, step_info; do_copy = true)
    if do_copy
        case = deepcopy(case)
    end
    (; model, forces, parameters, dt) = case
    is_multimodel = case.model isa MultiModel

    dd_dict = opt_dict[:model]
    domain = reservoir_domain(model)
    changed_dd = false
    for (k, v) in pairs(dd_dict)
        if haskey(domain, k)
            changed_dd = true
            domain[k] = v
        end
    end
    if changed_dd
        # Now we call the setup again
        parameters = setup_parameters(model)
    end
    if is_multimodel
        rparameters = parameters[:Reservoir]
    else
        rparameters = parameters
    end
    if haskey(dd_dict, :pore_volume)
        if haskey(rparameters, :StaticFluidVolume)
            k = :StaticFluidVolume
        else
            k = :FluidVolume
        end
        rparameters[k] = dd_dict[:pore_volume]
    end
    for (k, v) in pairs(opt_dict[:parameters])
        rparameters[k] = v
    end
    if is_multimodel
        for (w, wsub) in pairs(opt_dict[:wells])
            for (k, v) in pairs(wsub)
                @assert haskey(parameters[w], k)
                @assert size(parameters[w][k]) == size(v)
                parameters[w][k] = v
            end
        end
    end
    init = opt_dict[:state0]
    rmodel = reservoir_model(model)
    for (k, var) in pairs(rmodel.primary_variables)
        if var isa Jutul.FractionVariables
            val = init[k]
            for i in axes(val, 2)
                v = zero(eltype(val))
                for j in axes(val, 1)
                    v += val[j, i]
                end
                if v == 0.0
                    val[:, i] .= 1.0 / size(val, 1)
                else
                    val[:, i] ./= v
                end
            end
        end
    end
    state0 = setup_reservoir_state(model, init)
    new_case = JutulCase(model, dt, forces, parameters = parameters, state0 = state0)
    return new_case
end

function optimize_reservoir(dopt, objective, setup_fn = dopt.setup_function;
        simulator_arg = (output_substates = true, ),
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
        simulator_arg = (output_substates = true, ),
        kwarg...
    )
    if ismissing(setup_fn)
        error("Setup function was not found in DictParameters struct or as last positional argument.")
    end
    case0 = setup_fn(dopt.parameters, missing)
    sim, cfg = setup_reservoir_simulator(case0; info_level = -1, simulator_arg...)
    return Jutul.parameters_gradient(dopt, objective, setup_fn; simulator = sim, config = cfg, kwarg...)
end
