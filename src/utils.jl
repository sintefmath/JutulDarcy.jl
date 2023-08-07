export reservoir_model
reservoir_model(model) = model
reservoir_model(model::MultiModel) = model.models.Reservoir

reservoir_storage(model, storage) = storage
reservoir_storage(model::MultiModel, storage) = storage.Reservoir


export reservoir_domain
"""
    reservoir_domain(g; permeability = 9.869232667160130e-14, porosity = 0.1, kwarg...)

Set up a `DataDomain` instance for given mesh or other representation `g`.
`permeability` and `porosity` are then added to the domain. If scalars are
passed, they are expanded to cover all cells. Arrays are asserted to match all
cells. Permeability is either one value per cell (diagonal scalar), one value
per dimension given in each row (for a diagonal tensor) or a vector that
represents a compact full tensor representation (6 elements in 3D, 3 in 2D).
"""
function reservoir_domain(g; permeability = 9.869232667160130e-14, porosity = 0.1, kwarg...)
    return DataDomain(g; permeability = permeability, porosity = porosity, kwarg...)
end

export setup_reservoir_model
"""
    setup_reservoir_model(reservoir, system; wells = [], <keyword arguments>)
    setup_reservoir_model(reservoir, system; wells = [], context = DefaultContext(), reservoir_context = nothing, backend = :csc, <keyword arguments>)

Set up a reservoir `MultiModel` for a given reservoir `DataDomain` typically set
up from  [`reservoir_domain`](@ref) and an optional vector of wells that are
created using [`setup_vertical_well`](@ref) and  [`setup_well`](@ref).

The routine automatically sets up a facility and couples the wells with the
reservoir and that facility.
"""
function setup_reservoir_model(reservoir::DataDomain, system;
    wells = [],
    context = DefaultContext(),
    reservoir_context = nothing,
    general_ad = false,
    backend = :csc,
    split_wells = false,
    assemble_wells_together = true,
    parameters = Dict{Symbol, Any}(),
    kwarg...
    )
    # List of models (order matters)
    models = OrderedDict{Symbol, Jutul.AbstractSimulationModel}()
    reservoir_context, context = Jutul.select_contexts(
        backend; 
        main_context = reservoir_context,
        context = context,
        kwarg...
    )
    # We first set up the reservoir
    models[:Reservoir] = SimulationModel(
        reservoir,
        system,
        context = reservoir_context,
        general_ad = general_ad
    )
    # Then we set up all the wells
    mode = PredictionMode()
    if length(wells) > 0
        for w in wells
            if w isa SimpleWell
                well_context = reservoir_context
            else
                well_context = context
            end
            w_domain = DataDomain(w)
            wname = w.name
            models[wname] = SimulationModel(w_domain, system, context = well_context)
            if split_wells
                wg = WellGroup([wname])
                F = SimulationModel(wg, mode, context = context, data_domain = DataDomain(wg))
                models[Symbol(string(wname)*string(:_ctrl))] = F
            end
        end
        # Add facility that gorups the wells
        if !split_wells
            wg = WellGroup(map(x -> x.name, wells))
            F = SimulationModel(wg, mode, context = context, data_domain = DataDomain(wg))
            models[:Facility] = F
        end
    end

    # Put it all together as multimodel
    model = reservoir_multimodel(models, split_wells = split_wells, assemble_wells_together = assemble_wells_together)
    # Insert domain here.
    parameters = setup_parameters(model, parameters)
    return (model, parameters)
end

export setup_reservoir_simulator
"""
    setup_reservoir_simulator(models, initializer, parameters = nothing; <keyword arguments>)

# Arguments
- `models`: either a single model or a Dict with the key :Reservoir for multimodels
- `initializer`: used to setup state0, must be compatible with `model`
- `parameters`: initialized parameters, must be compatible with `model` if provided
- `linear_solver=:bicgstab`: iterative solver to use (provided model supports it)
- `precond=:cpr`: preconditioner for iterative solver: Either :cpr or :ilu0.
- `rtol=1e-3`: relative tolerance for linear solver
- `initial_dt=3600*24.0`: initial time-step in seconds (one day by default)
- `target_its=8`: target number of nonlinear iterations per time step
- `offset_its=1`: dampening parameter for time step selector where larger values lead to more pessimistic estimates.
- `tol_cnv=1e-3`: maximum allowable point-wise error (volume-balance)
- `tol_mb=1e-7`: maximum alllowable integrated error (mass-balance)
- `specialize=false`: use deep specialization of storage for faster execution, but significantly more compile time

Additional keyword arguments are passed onto [`simulator_config`](@ref).
"""
function setup_reservoir_simulator(models, initializer, parameters = nothing;
                                                        specialize = false,
                                                        split_wells = false,
                                                        assemble_wells_together = true,
                                                        kwarg...)
    if isa(models, SimulationModel)
        DT = Dict{Symbol, Any}
        models = DT(:Reservoir => models)
        initializer = DT(:Reservoir => initializer)
        if !isnothing(parameters)
            parameters = DT(:Reservoir => parameters)
        end
    end
    # Convert to multi model
    mmodel = reservoir_multimodel(models, specialize = specialize, split_wells = split_wells, assemble_wells_together = assemble_wells_together)
    if isnothing(parameters)
        parameters = setup_parameters(mmodel)
    end
    # Set up simulator itself, containing the initial state
    state0 = setup_state(mmodel, initializer)

    case = JutulCase(mmodel, state0 = state0, parameters = parameters)
    setup_reservoir_simulator(case; kwarg...)
end

function mode_to_backend(mode::Symbol)
    if mode == :mpi
        mode = MPI_PArrayBackend()
    elseif mode == :parray
        mode = Jutul.JuliaPArrayBackend()
    else
        @assert mode == :debug "Mode must be one of :mpi, :parray, :debug, was :$mode"
        mode = Jutul.DebugPArrayBackend()
    end
end

function mode_to_backend(mode::Jutul.PArrayBackend)
    return mode
end

function setup_reservoir_simulator(case::JutulCase;
                            mode = :default,
                            precond = :cpr,
                            linear_solver = :bicgstab,
                            max_dt = Inf,
                            rtol = nothing,
                            initial_dt = 3600.0*24.0,
                            target_ds = Inf,
                            target_its = 8,
                            offset_its = 1,
                            tol_cnv = 1e-3,
                            tol_mb = 1e-7,
                            info_level = 0,
                            tol_cnv_well = 10*tol_cnv,
                            tol_mb_well = 1e4*tol_mb,
                            set_linear_solver = linear_solver isa Symbol,
                            timesteps = :auto,
                            parray_arg = NamedTuple(),
                            linear_solver_arg = NamedTuple(),
                            extra_timing_setup = false,
                            kwarg...)
    if mode == :default
        sim = Simulator(case, extra_timing = extra_timing_setup)
    else
        b = mode_to_backend(mode)
        sim = setup_reservoir_simulator_parray(case, b; parray_arg...);
    end
    t_base = TimestepSelector(initial_absolute = initial_dt, max = max_dt)
    sel = Vector{Any}()
    push!(sel, t_base)
    if timesteps == :auto || timesteps == :iteration
        if isfinite(target_its)
            t_its = IterationTimestepSelector(target_its, offset = offset_its)
            push!(sel, t_its)
        end

        if isfinite(target_ds)
            @assert mode == :default "target_ds is only supported in serial."
            t_sat = VariableChangeTimestepSelector(
                :Saturations, target_ds, relative = false, reduction = :max, model = :Reservoir
                )
            push!(sel, t_sat)
        end
    else
        @assert isnothing(timesteps)
        push!(sel, t_base)
    end
    # Config: Linear solver, timestep selection defaults, etc...
    if set_linear_solver
        if info_level < 1
            v = -1
        else
            v = 0
        end
        lsolve = reservoir_linsolve(case.model, precond;
                                            rtol = rtol,
                                            linear_solver_arg...,
                                            )
        extra_arg = (linear_solver = lsolve, )
    else
        extra_arg = NamedTuple()
    end
    cfg = simulator_config(sim; extra_arg..., timestep_selectors = sel, info_level = info_level, kwarg...)
    set_default_cnv_mb!(cfg, sim, tol_cnv = tol_cnv, tol_mb = tol_mb, tol_cnv_well = tol_cnv_well, tol_mb_well = tol_mb_well)
    return (sim, cfg)
end

export simulate_reservoir

function simulate_reservoir(state0, model, dt;
        parameters = setup_parameters(model),
        restart = false,
        forces = setup_forces(model),
        kwarg...
    )
    sim, config = setup_reservoir_simulator(model, state0, parameters; kwarg...)
    result = simulate!(sim, dt, forces = forces, config = config, restart = restart);
    return ReservoirSimResult(model, result, forces)
end

function simulate_reservoir(case::JutulCase; restart = false, kwarg...)
    (; model, forces, state0, parameters, dt) = case
    sim, config = setup_reservoir_simulator(model, state0, parameters; kwarg...)
    result = simulate!(sim, dt, forces = forces, config = config, restart = restart);
    return ReservoirSimResult(model, result, forces)
end

function set_default_cnv_mb!(cfg::JutulConfig, sim::JutulSimulator; kwarg...)
    set_default_cnv_mb!(cfg, sim.model; kwarg...)
end

function set_default_cnv_mb!(cfg, model; kwarg...)
    set_default_cnv_mb_inner!(cfg[:tolerances], model; kwarg...)
end

function set_default_cnv_mb_inner!(tol, model; tol_cnv = 1e-3, tol_mb = 1e-7, tol_mb_well = 1e-3, tol_cnv_well = 1e-2)
    sys = model.system
    if sys isa ImmiscibleSystem || sys isa BlackOilSystem || sys isa CompositionalSystem
        if physical_representation(model) isa WellDomain
            c = tol_cnv_well
            m = tol_mb_well
        else
            c = tol_cnv
            m = tol_mb
        end
        tol[:mass_conservation] = (CNV = c, MB = m)
    end
end

function set_default_cnv_mb!(cfg, model::MultiModel; kwarg...)
    for (k, m) in pairs(model.models)
        set_default_cnv_mb_inner!(cfg[:tolerances][k], m; kwarg...)
    end
end

function setup_reservoir_cross_terms!(model::MultiModel)
    rmodel = reservoir_model(model)
    has_composite = rmodel isa Jutul.CompositeModel
    if has_composite
        systems = rmodel.system.systems
        has_flow = haskey(systems, :flow)
        has_thermal = haskey(systems, :thermal)
        conservation = Pair(:flow, :mass_conservation)
        energy = Pair(:thermal, :energy_conservation)
    else
        has_flow = rmodel.system isa MultiPhaseSystem
        has_thermal = !has_flow
        conservation = :mass_conservation
        energy = :energy_conservation
    end
    for (k, m) in pairs(model.models)
        if k == :Reservoir
            # These are set up from wells via symmetry
        elseif m.domain isa WellGroup
            for target_well in m.domain.well_symbols
                if has_flow
                    ct = FacilityFromWellFlowCT(target_well)
                    add_cross_term!(model, ct, target = k, source = target_well, equation = :control_equation)

                    ct = WellFromFacilityFlowCT(target_well)
                    add_cross_term!(model, ct, target = target_well, source = k, equation = conservation)
                end
                if has_thermal
                    ct = WellFromFacilityThermalCT(target_well)
                    add_cross_term!(model, ct, target = target_well, source = k, equation = energy)
                end
            end
        else
            g = physical_representation(m.domain)
            if g isa WellDomain
                WI = vec(g.perforations.WI)
                rc = vec(g.perforations.reservoir)
                wc = vec(g.perforations.self)
                # Put these over in cross term
                if has_flow
                    ct = ReservoirFromWellFlowCT(WI, rc, wc)
                    add_cross_term!(model, ct, target = :Reservoir, source = k, equation = conservation)
                end
                if has_thermal
                    CI = 1000 .* WI
                    ct = ReservoirFromWellThermalCT(CI, WI, rc, wc)
                    add_cross_term!(model, ct, target = :Reservoir, source = k, equation = energy)
                end
            end
        end
    end
end

function reservoir_multimodel(model::MultiModel; kwarg...)
    # The multimodel is a reservoir multimodel if there exists a submodel named Reservoir
    @assert haskey(model.models, :Reservoir)
    return model
end

function reservoir_multimodel(models::AbstractDict; specialize = false, split_wells = false, assemble_wells_together = haskey(models, :Facility))
    res_model = models[:Reservoir]
    is_block(x) = Jutul.is_cell_major(matrix_layout(x.context))
    block_backend = is_block(res_model)
    n = length(models)
    if block_backend && n > 1
        if  !(split_wells == true) || assemble_wells_together
            groups = repeat([2], n)
            for (i, k) in enumerate(keys(models))
                m = models[k]
                if is_block(m)
                    groups[i] = 1
                end
            end
        else
            groups = repeat([1], n)
            gpos = 1
            pos = 2
            wno = 1
            mkeys = collect(keys(models))
            nw = 1
            for (k, m) in models
                if endswith(String(k), "_ctrl")
                    nw += 1
                end
            end
            if split_wells == :threads
                npartition = Threads.nthreads()
            else
                npartition = nw
            end
            p = Jutul.partition_linear(npartition, nw)
            for (k, m) in models
                if k != :Reservoir
                    sk = String(k)
                    if endswith(sk, "_ctrl")
                        wk = Symbol(sk[1:end-5])
                        wpos = findfirst(isequal(wk), mkeys)
                        if false
                            g = pos
                        else
                            g = p[wno]+1
                        end
                        groups[wpos] = g
                        groups[gpos] = g
                        pos += 1
                        wno += 1
                    end
                end
                gpos += 1
            end
        end
        red = :schur_apply
        outer_context = DefaultContext()
    else
        outer_context = models[:Reservoir].context
        groups = nothing
        red = nothing
    end
    models = convert_to_immutable_storage(models)
    model = MultiModel(models, groups = groups, context = outer_context, reduction = red, specialize = specialize)
    setup_reservoir_cross_terms!(model)
    return model
end

export setup_reservoir_state
"""
    setup_reservoir_state(model, <keyword arguments>)
    # Ex: For immiscible two-phase
    setup_reservoir_state(model, Pressure = 1e5, Saturations = [0.2, 0.8])

Convenience constructor that initializes a state for a `MultiModel` set up using
[`setup_reservoir_model`](@ref). The main convenience over [`setup_state`](@ref)
is only the reservoir initialization values need be provided: wells are
automatically initialized from the connected reservoir cells.

As an alternative to passing keyword arguments, a `Dict{Symbol, Any}` instance
can be sent in as a second, non-keyword argument.
"""
function setup_reservoir_state(model::MultiModel; kwarg...)
    rmodel = reservoir_model(model)
    pvars = collect(keys(Jutul.get_primary_variables(rmodel)))
    res_state = setup_reservoir_state(rmodel; kwarg...)
    # Next, we initialize the wells.
    init = Dict(:Reservoir => res_state)
    perf_subset(v::AbstractVector, i) = v[i]
    perf_subset(v::AbstractMatrix, i) = v[:, i]
    perf_subset(v, i) = v
    for k in keys(model.models)
        if k == :Reservoir
            # Already done
            continue
        end
        W = model.models[k]
        if W.domain isa WellGroup
            # Facility or well group
            init_w = setup_state(W, TotalSurfaceMassRate = 0.0)
        else
            # Wells
            init_w = Dict{Symbol, Any}()
            W = model.models[k]
            wg = physical_representation(W.domain)
            res_c = wg.perforations.reservoir
            if wg isa MultiSegmentWell
                # Repeat top node. Not fully robust.
                c = res_c[vcat(1, 1:length(res_c))]
                init_w[:TotalMassFlux] = 0.0
            else
                c = res_c[1]
            end
            for pk in pvars
                pv = res_state[pk]
                init_w[pk] = perf_subset(pv, c)
            end
        end
        init[k] = init_w
    end
    state = setup_state(model, init)
    return state
end

function setup_reservoir_state(model, init)
    return setup_reservoir_state(model; pairs(init)...)
end

function setup_reservoir_state(rmodel::SimulationModel; kwarg...)
    pvars = collect(keys(Jutul.get_primary_variables(rmodel)))
    svars = collect(keys(Jutul.get_secondary_variables(rmodel)))
    np = length(pvars)
    found = Symbol[]
    res_init = Dict{Symbol, Any}()
    for (k, v) in kwarg
        I = findfirst(isequal(k), pvars)
        if isnothing(I)
            if !(k in svars)
                @warn "Recieved primary variable $k, but this is not known to reservoir model... Adding anyway."
            end
        else
            push!(found, k)
        end
        res_init[k] = v
    end
    handle_alternate_primary_variable_spec!(res_init, found, rmodel.system)
    if length(found) != length(pvars)
        missing_primary_variables = setdiff(pvars, found)
        @warn "Not all primary variables were initialized for reservoir model." missing_primary_variables
    end
    return setup_state(rmodel, res_init)
end

function handle_alternate_primary_variable_spec!(res_init, found, system)
    # Internal utility to handle non-trivial specification of primary variables
    return res_init
end

export setup_reservoir_forces
"""
    setup_reservoir_forces(model; control = nothing, limits = nothing, set_default_limits = true, <keyword arguments>)

Set up driving forces for a reservoir model with wells
"""
function setup_reservoir_forces(model::MultiModel; control = nothing, limits = nothing, set_default_limits = true, kwarg...)
    submodels = model.models
    has_facility = any(x -> isa(x.domain, WellGroup), values(submodels))
    no_well_controls = isnothing(control) && isnothing(limits)
    @assert no_well_controls || has_facility "Model must have facility."
    if haskey(submodels, :Facility)
        # Unified facility for all wells
        facility = model.models.Facility

        surface_forces = setup_forces(facility,
            control = control,
            limits = limits,
            set_default_limits = set_default_limits
        )
        # Set up forces for the whole model.
        return setup_forces(model, Facility = surface_forces; kwarg...)
    else
        new_forces = Dict{Symbol, Any}()
        for (k, m) in pairs(submodels)
            if m isa SimpleWellFlowModel || m isa MSWellFlowModel
                ctrl_symbol = Symbol("$(k)_ctrl")
                @assert haskey(submodels, ctrl_symbol) "Controller for well $k must be present with the name $ctrl_symbol"
                facility = submodels[ctrl_symbol]
                new_forces[k] = setup_forces(facility,
                    control = control,
                    limits = limits,
                    set_default_limits = set_default_limits
                )
            end
        end
        return setup_forces(model; pairs(new_forces)..., kwarg...)
    end
end

export full_well_outputs, well_output, well_symbols, wellgroup_symbols, available_well_targets

"""
    full_well_outputs(model, states, forces; targets = available_well_targets(model.models.Reservoir), shortname = false)

Get the full set of well outputs after a simulation has occured, for plotting or other post-processing.
"""
function full_well_outputs(model, states, forces; targets = available_well_targets(model.models.Reservoir), shortname = false)
    out = Dict()
    if shortname
        tm = :mass
    else
        tm = Symbol("Total surface mass rate")
    end
    for w in well_symbols(model)
        out[w] = Dict()
        for t in targets
            out[w][translate_target_to_symbol(t(1.0), shortname = shortname)] = well_output(model, states, w, forces, t)
        end
        out[w][Symbol(tm)] = well_output(model, states, w, forces, :TotalSurfaceMassRate)
    end
    return out
end

"""
    well_output(model, states, well_symbol, forces, target = BottomHolePressureTarget)

Get a specific well output from a valid operational target once a simulation is completed an `states` are available.
"""
function well_output(model::MultiModel, states, well_symbol, forces, target = BottomHolePressureTarget)
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
    rhoS_o = reference_densities(model.models[well_symbol].system)

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
        if forces isa AbstractVector
            force = forces[i]
        else
            force = forces
        end
        if haskey(force, :outer)
            force = force.outer
        end
        if target == :TotalSurfaceMassRate
            d[i] = q_t
        else
            if q_t == 0
                current_control = DisabledControl()
                if target == BottomHolePressureTarget
                    v = well_state.Pressure[1]
                else
                    v = 0.0
                end
                d[i] = v
            else
                if haskey(force, :Facility)
                    gforce = force[:Facility]
                else
                    gforce = force[Symbol("$(well_symbol)_ctrl")]
                end
                control = gforce.control[well_symbol]
                current_control = replace_target(control, BottomHolePressureTarget(1.0))
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
        if isa(physical_representation(D), WellDomain)
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

function partitioner_input(model, parameters; conn = :trans)
    rmodel = reservoir_model(model)
    if haskey(parameters, :Reservoir)
        parameters = parameters[:Reservoir]
    end
    grid = physical_representation(rmodel.domain)

    N = grid.neighborship
    trans = parameters[:Transmissibilities]
    if conn == :trans
        T = copy(trans)
    else
        @assert conn == :unit
        T = ones(Int, length(trans))
    end
    groups = []
    if model isa MultiModel
        for (k, m) in pairs(model.models)
            wg = physical_representation(m.domain)
            if wg isa WellDomain
                push!(groups, copy(wg.perforations.reservoir))
            end
        end
    end
    return (N, T, groups)
end

function reservoir_partition(model::MultiModel, p)
    p_res = SimplePartition(p)
    models = model.models
    function model_is_well(m)
        d = physical_representation(m.domain)
        return isa(d, JutulDarcy.WellDomain)
    end
    wpart = Dict()
    for key in keys(models)
        m = models[key]
        if model_is_well(m)
            wg = physical_representation(m.domain)
            wc = wg.perforations.reservoir
            unique_part_wcells = unique(p[wc])
            @assert length(unique_part_wcells) == 1 "All cells of well $key must be in the same partitioned block. Found: $unique_part_wcells for cells $wc = $(p[wc])"
            wpart[key] = first(unique_part_wcells)
        end
    end
    part = Dict{Symbol, Any}()
    for key in keys(models)
        m = models[key]
        if key == :Reservoir
            part[key] = p_res
        elseif model_is_well(m)
            part[key] = wpart[key]
        else
            # Well group, probably?
            s = Symbol(string(key)[1:end-5])
            part[key] = wpart[s]
        end
    end
    return SimpleMultiModelPartition(part, :Reservoir)
end

function reservoir_partition(model, p)
    return SimplePartition(p)
end

function Base.iterate(t::ReservoirSimResult)
    return (t.wells, :wells)
end

function Base.iterate(t::ReservoirSimResult, state)
    if state == :states
        return (t.time, nothing)
    else
        @assert state == :wells "Unapck syntax: ws, states, t = res_result"
        return (t.states, :states)
    end
end

Base.lastindex(t::ReservoirSimResult) = :states

function Base.show(io::IO, ::MIME"text/plain", sr::ReservoirSimResult)
    # return 
    function print_keys(prefix, el)
        for k in keys(el)
            v = el[k]
            if v isa AbstractDict
                print(io, "$prefix:$k\n")
                print_keys("  $prefix", v)
            else
                if v isa AbstractVecOrMat
                    s = " of size $(size(v))"
                else
                    s = ""
                end
                print(io, "$prefix:$k => $(typeof(v))$s\n")
            end
        end
    end
    # fmt = raw"u. dd Y H:mm"

    states = sr.states
    n = length(states)
    print(io, sr)
    print(io, ":\n")
    if n > 0
        wk = keys(sr.wells)
        nw = length(wk)
        print(io, "\n  wells ($nw present):\n")
        if nw > 0
            for k in wk
                print(io, "    :$k\n")
            end
            print(io, "    Results per well:\n")
            print_keys("       ", sr.wells[first(wk)])
        end
        el = first(states)
        print(io, "\n  states (Vector with $n entries, reservoir variables for each state)\n")
        print_keys("    ", el)
    end
    print(io, "\n  time (report time for each state)\n     $(typeof(sr.time)) of length $n\n")
    print(io, "\n  result (extended states, reports)\n     $(sr.result)\n")
    extra_keys = keys(sr.extra)
    print(io, "\n  extra\n     $(typeof(sr.extra)) with ")
    if length(extra_keys) == 0
        print(io, "with no data.\n")
    else
        ek = join(extra_keys, ", :")
        print(io, "keys :$ek\n")
    end
    Jutul.print_sim_result_timing(io, sr.result)
end

function Base.show(io::IO, sr::ReservoirSimResult)
    n = length(sr.states)
    if n == 1
        s = "entry"
    else
        s = "entries"
    end
    print(io, "ReservoirSimResult with $n $s")
end

"""
    generate_phase_indices(phases, canonical = (a = AqueousPhase(), l = LiquidPhase(), v = VaporPhase()))

Generate mapping for canonical ordering of phases.
"""
function generate_phase_indices(phases, canonical = (a = AqueousPhase(), l = LiquidPhase(), v = VaporPhase()))
    phase_ind_dict = OrderedDict()
    for (s, ph) in pairs(canonical)
        pos = findfirst(isequal(ph), phases)
        if !isnothing(pos)
            phase_ind_dict[s] = pos
        end
    end
    return NamedTuple(pairs(phase_ind_dict))
end