function setup_reservoir_dict_optimization(dict, F = missing; kwarg...)
    return DictParameters(dict, F; kwarg...)
end

function setup_reservoir_dict_optimization(case::JutulCase;
        use_trans = false,
        use_pore_volume = false,
        use_multipliers = false,
        strict = false,
        verbose = true,
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
    ncells = number_of_cells(rdomain)
    nfaces = number_of_faces(rdomain)

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
        trans = copy(rparameters[:Transmissibilities])
        if use_multipliers
            prm_dict[:multiplier_trans] = ones(nfaces)
        else
            trans::Vector{Float64}
            for (i, v) in enumerate(trans)
                # Avoid for optimizer
                trans[i] = max(v, 1e-20)
            end
            prm_dict[:Transmissibilities] = trans
        end
        push!(skip_list, :permeability)
    elseif use_multipliers
        prm_dict[:multiplier_permeability] = ones(ncells)
        prm_dict[:multiplier_vertical_permeability] = ones(ncells)
        push!(skip_list, :permeability)
    end
    if use_pore_volume || use_multipliers
        if haskey(rparameters, :StaticFluidVolume)
            k = :StaticFluidVolume
        else
            k = :FluidVolume
        end
        if use_multipliers
            prm_dict[:multiplier_pore_volume] = ones(ncells)
        else
            dd_dict[:pore_volume] = copy(rparameters[k])
        end
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
                wi = copy(subprm[:WellIndices])
                if use_multipliers
                    subdict[:multiplier_wellindices] = ones(length(wi))
                else
                    subdict[:WellIndices] = copy(subprm[:WellIndices])
                end
                if haskey(subprm, :WellIndicesThermal) && is_thermal
                    if use_multipliers
                        subdict[:multiplier_wellindices_thermal] = ones(length(wi))
                    else
                        subdict[:WellIndicesThermal] = copy(subprm[:WellIndicesThermal])
                    end
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
    F(D, step_info = missing) = optimization_resetup_reservoir_case(D, case, step_info,
        do_copy = do_copy,
        use_multipliers = use_multipliers,
        use_pore_volume = use_pore_volume,
        use_trans = use_trans
    )
    return DictParameters(opt_dict, F, strict = strict, verbose = verbose)
end

function optimization_resetup_reservoir_case(opt_dict::AbstractDict, case::JutulCase, step_info;
        do_copy = true,
        use_trans = false,
        use_pore_volume = false,
        use_multipliers = false
    )
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
            e = associated_entity(domain, k)
            changed_dd = true
            domain[k, e] = v
        end
    end
    if use_multipliers
        if !use_trans
            perm = domain[:permeability]
            permmult = opt_dict[:parameters][:multiplier_permeability]
            perm = perm .* permmult'
            if size(perm, 1) > 2
                vpermmult = opt_dict[:parameters][:multiplier_vertical_permeability]
                for i in axes(perm, 2)
                    perm[3, i] *= vpermmult[i]
                end
            end
            domain[:permeability] = perm
            changed_dd = true
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
    if use_pore_volume || use_multipliers
        if haskey(rparameters, :StaticFluidVolume)
            k = :StaticFluidVolume
        else
            k = :FluidVolume
        end
        if use_multipliers
            pv_val = rparameters[k] .* opt_dict[:parameters][:multiplier_pore_volume]
        else
            pv_val = dd_dict[:pore_volume]
        end
        rparameters[k] = pv_val
    end
    for (k, v) in pairs(opt_dict[:parameters])
        rparameters[k] = v
    end
    if is_multimodel
        for (w, wsub) in pairs(opt_dict[:wells])
            wprm = parameters[w]
            for (k, v) in pairs(wsub)
                if k == :multiplier_wellindices
                    wi = wprm[:WellIndices]
                    wprm[:WellIndices] = wi .* v
                elseif k == :multiplier_wellindices_thermal
                    wi = wprm[:WellIndicesThermal]
                    wprm[:WellIndicesThermal] = wi .* v
                else
                    @assert haskey(wprm, k)
                    @assert size(wprm[k]) == size(v)
                    wprm[k] = v
                end
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

"""
    optimize_reservoir(dopt, objective, setup_fn = dopt.setup_function)

Perform optimization for reservoir models using the given `DictParameters`
struct `dopt`, objective function `objective`, and setup function `setup_fn`.

Additional keyword arguments can be provided to customize the simulator
setup and optimization process:
- `simulator_arg`: Arguments for the simulator setup (default: `(output_substates = true,)`).
- `simulator`: Custom simulator instance (default: `missing`).
- `config`: Configuration for the simulator (default: `missing`).
- `deps`: Dependencies for the optimization (default: `:parameters_and_state0`).

# Notes
If you are optimizing forces (i.e. well constraints or boundary conditions), you
need to set `deps = :case` to ensure that gradients are correctly computed. The
same applies if you change the model itself in the setup function. The defaults
of `:parameters_and_state0` provide significant performance improvements in most cases
where only parameters and initial state are changed.
"""
function optimize_reservoir(dopt, objective, setup_fn = dopt.setup_function;
        info_level = 0,
        simulator_arg = (output_substates = true, info_level = info_level, end_report = info_level > 0),
        simulator = missing,
        config = missing,
        deps = :parameters_and_state0,
        kwarg...
    )
    sim, cfg = setup_simulator_for_reservoir_optimization(dopt, setup_fn, simulator, config, simulator_arg)
    return Jutul.optimize(dopt, objective, setup_fn; simulator = sim, config = cfg, deps = deps, info_level = info_level, kwarg...)
end

function parameters_gradient_reservoir(dopt, objective, setup_fn = dopt.setup_function;
        simulator_arg = (output_substates = true, ),
        simulator = missing,
        config = missing,
        kwarg...
    )
    sim, cfg = setup_simulator_for_reservoir_optimization(dopt, setup_fn, simulator, config, simulator_arg)
    return Jutul.parameters_gradient(dopt, objective, setup_fn; simulator = sim, config = cfg, kwarg...)
end

function setup_simulator_for_reservoir_optimization(dopt, setup_fn, simulator, config, simulator_arg)
    if ismissing(setup_fn)
        error("Setup function was not found in DictParameters struct or as last positional argument.")
    end
    has_sim = !ismissing(simulator)
    has_cfg = !ismissing(config)
    if !has_sim && !has_cfg
        case0 = setup_fn(dopt.parameters, missing)
        simulator, config = setup_reservoir_simulator(case0; info_level = -1, simulator_arg...)
    else
        has_sim == has_cfg || error("Simulator and config must be provided together")
    end
    return (simulator, config)
end
