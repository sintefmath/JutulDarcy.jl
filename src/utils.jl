reservoir_model(model) = model
reservoir_model(model::MultiModel) = model.models.Reservoir

reservoir_storage(model, storage) = storage
reservoir_storage(model::MultiModel, storage) = storage.Reservoir


export setup_reservoir_simulator
function setup_reservoir_simulator(models, initializer, parameters = nothing; method = :cpr, rtol = 0.005, initial_dt = 3600.0*24.0, target_its = 8, offset_its = 1, kwarg...)
    res_model = models[:Reservoir]
    # Convert to multi model
    block_backend = is_cell_major(matrix_layout(res_model.context))
    if block_backend && length(models) > 1
        groups = repeat([2], length(models))
        groups[1] = 1
        red = :schur_apply
        outer_context = DefaultContext()
    else
        outer_context = models[:Reservoir].context
        groups = nothing
        red = nothing
    end
    mmodel = MultiModel(convert_to_immutable_storage(models), groups = groups, 
                                                              context = outer_context,
                                                              reduction = red)
    # Set up simulator itself, containing the initial state
    state0 = setup_state(mmodel, initializer)
    sim = Simulator(mmodel, state0 = state0, parameters = deepcopy(parameters))

    # Config: Linear solver, timestep selection defaults, etc...
    lsolve = reservoir_linsolve(mmodel, method, rtol = rtol)
    # day = 3600.0*24.0
    t_base = TimestepSelector(initial_absolute = initial_dt, max = Inf)
    t_its = IterationTimestepSelector(target_its, offset = offset_its)
    cfg = simulator_config(sim, timestep_selectors = [t_base, t_its], linear_solver = lsolve; kwarg...)

    return (sim, cfg)
end

export full_well_outputs, well_output, well_symbols, wellgroup_symbols, available_well_targets

function full_well_outputs(model, parameters, states; targets = available_well_targets(model.models.Reservoir), shortname = false)
    out = Dict()
    if shortname
        tm = "MASS"
    else
        tm = "Total surface mass rate"
    end
    for w in well_symbols(model)
        out[w] = Dict()
        for t in targets
            out[w][translate_target_to_symbol(t(1.0), shortname = shortname)] = well_output(model, parameters, states, w, t)
        end
        out[w][Symbol(tm)] = well_output(model, parameters, states, w, :TotalSurfaceMassRate)
    end
    return out
end

function well_output(model::MultiModel, parameters, states, well_symbol, target = BottomHolePressureTarget)
    n = length(states)
    d = zeros(n)

    groups = wellgroup_symbols(model)
    group = nothing
    for g in groups
        if well_symbol in model.models[g].domain.well_symbols
            group = g
            break
        end
    end
    rhoS_o = parameters[well_symbol][:reference_densities]

    to_target(t::DataType) = t(1.0)
    to_target(t::Type) = t(1.0)
    to_target(t::Symbol) = t

    target_limit = to_target(target)

    pos = get_well_position(model.models[group].domain, well_symbol)
    well_model = model.models[well_symbol]
    for (i, state) = enumerate(states)
        well_state = state[well_symbol]
        well_state = convert_to_immutable_storage(well_state)
        q_t = state[group][:TotalSurfaceMassRate][pos]
        if target == :TotalSurfaceMassRate
            d[i] = q_t
        else
            if q_t == 0
                current_control = DisabledControl()
                d[i] = 0.0
            else
                if q_t < 0
                    current_control = ProducerControl(BottomHolePressureTarget(1.0))
                else
                    current_control = InjectorControl(BottomHolePressureTarget(1.0), 1.0)
                end
                rhoS, S = flash_wellstream_at_surface(well_model, well_state, rhoS_o)
                v = well_target_value(q_t, current_control, target_limit, well_model, well_state, rhoS, S)
                d[i] = v
            end
        end
    end
    return d
end

function well_symbols(model::MultiModel)
    models = model.models
    symbols = Vector{Symbol}()
    for (k, m) in pairs(models)
        D = m.domain
        if isa(D, DiscretizedDomain) && isa(D.grid, WellGrid)
            push!(symbols, k)
        end
    end
    return symbols
end

function wellgroup_symbols(model::MultiModel)
    models = model.models
    symbols = Vector{Symbol}()
    for (k, m) in pairs(models)
        D = m.domain
        if isa(D, WellControllerDomain)
            push!(symbols, k)
        end
    end
    return symbols
end

function available_well_targets(model)
    phases = get_phases(model.system)
    targets = [BottomHolePressureTarget, SurfaceLiquidRateTarget, TotalRateTarget]
    if AqueousPhase() in phases
        push!(targets, SurfaceLiquidRateTarget)
        push!(targets, SurfaceWaterRateTarget)
    end
    if LiquidPhase() in phases
        push!(targets, SurfaceLiquidRateTarget)
        push!(targets, SurfaceOilRateTarget)
    end
    if VaporPhase() in phases
        push!(targets, SurfaceGasRateTarget)
    end
    return unique(targets)
end
