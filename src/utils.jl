"""
    reservoir_model(model)

Get the reservoir model from a `MultiModel` or return the model itself if it is
not a `MultiModel`.
"""
function reservoir_model

end


function reservoir_model(model::SimulationModel; kwarg...)
    return model
end

function reservoir_model(model::MultiModel; kwarg...)
    return reservoir_model(model.models.Reservoir; kwarg...)
end

function reservoir_model(case::JutulCase; kwarg...)
    return reservoir_model(case.model; kwarg...)
end

function reservoir_model(model::Jutul.CompositeModel; type = missing)
    if !ismissing(type)
        model = Jutul.composite_submodel(model, type)
    end
    return model
end

reservoir_storage(model, storage) = storage
reservoir_storage(model::MultiModel, storage) = storage.Reservoir


"""
    reservoir_domain(g; permeability = convert_to_si(0.1, :darcy), porosity = 0.1, kwarg...)

Set up a `DataDomain` instance for given mesh or other representation `g`.
`permeability` and `porosity` are then added to the domain. If scalars are
passed, they are expanded to cover all cells. Arrays are asserted to match all
cells. Permeability is either one value per cell (diagonal scalar), one value
per dimension given in each row (for a diagonal tensor) or a vector that
represents a compact full tensor representation (6 elements in 3D, 3 in 2D).

# Default data and their values

| Name                         | Explanation                                | Unit         | Default |
|------------------------------|--------------------------------------------|--------------|---------|
| `permeability`               | Rock ability to conduct fluid flow         | ``m^2``      | 100 mD  |
| `porosity`                   | Rock void fraction open to flow (0 to 1)   | -            |  0.3    |
| `rock_density`               | Mass density of rock                       | ``kg^3/m^3`` | 2000.0  |
| `rock_heat_capacity`         | Specific heat capacity of rock             | ``J/(kg K)`` | 900.0   |
| `rock_thermal_conductivity`  | Heat conductivity of rock                  | ``W/m K``    | 3.0     |
| `fluid_thermal_conductivity` | Heat conductivity of fluid phases          | ``W/m K``    | 0.6     |
| `component_heat_capacity`    | Specific heat capacity of fluid components | ``J/(kg K)`` | 4184.0  |

Note that the default values are taken to be roughly those of water for fluid
phases and sandstone for those of rock. Choice of values can severely impact
your simulation results - take care to check the values that your physical
system makes use of!
"""
function reservoir_domain(g;
        permeability = convert_to_si(0.1, :darcy),
        porosity = 0.1,
        rock_thermal_conductivity = 3.0, # W/m K (~sandstone)
        fluid_thermal_conductivity = 0.6, # W/m K (~water)
        rock_heat_capacity = 900.0, # ~sandstone
        component_heat_capacity = 4184.0, # ~water
        rock_density = 2000.0,
        diffusion = missing,
        kwarg...
    )
    all(isfinite, permeability) || throw(ArgumentError("Keyword argument permeability has non-finite entries."))
    all(isfinite, porosity) || throw(ArgumentError("Keyword argument porosity has non-finite entries."))
    all(isfinite, fluid_thermal_conductivity) || throw(ArgumentError("Keyword argument fluid_thermal_conductivity has non-finite entries."))
    all(isfinite, rock_heat_capacity) || throw(ArgumentError("Keyword argument rock_heat_capacity has non-finite entries."))
    all(isfinite, component_heat_capacity) || throw(ArgumentError("Keyword argument component_heat_capacity has non-finite entries."))

    if !ismissing(diffusion)
        all(isfinite, diffusion) || throw(ArgumentError("Keyword argument diffusion has non-finite entries."))
        kwarg = (diffusion = diffusion, kwarg...)
    end
    minimum(permeability) >= 0 || throw(ArgumentError("All permeability values must be non-negative."))
    nk = length(permeability)
    nc = number_of_cells(g)
    if nk != nc && permeability isa AbstractVector
        d = dim(g)
        if nk == d || (d == 2 && nk == 3) || (d == 3 && nk == 6)
            permeability = repeat(permeability, 1, nc)
        end
    end

    reservoir = DataDomain(g;
        permeability = permeability,
        porosity = porosity,
        rock_thermal_conductivity = rock_thermal_conductivity,
        fluid_thermal_conductivity = fluid_thermal_conductivity,
        rock_density = rock_density,
        kwarg...
    )
    for k in [:porosity, :net_to_gross]
        if haskey(reservoir, k)
            val = reservoir[k]
            minimum(val) > 0 || throw(ArgumentError("Keyword argument $k must have positive entries."))
        end
    end
    return reservoir
end

"""
    reservoir_domain(m::Union{SimulationModel, MultiModel})

Get reservoir domain embedded in model.
"""
function reservoir_domain(m::Union{SimulationModel, MultiModel})
    d = reservoir_model(m).data_domain
    return d::DataDomain
end

"""
    reservoir_domain(case::JutulCase)

Get reservoir domain from a reservoir simulation case.
"""
function reservoir_domain(case::JutulCase)
    return reservoir_domain(case.model)
end

export reservoir_system

function reservoir_system(flow::MultiPhaseSystem; kwarg...)
    reservoir_system(;flow = flow, kwarg...)
end

export get_model_wells

function get_model_wells(case::JutulCase)
    return get_model_wells(case.model)
end

"""
    get_model_wells(model_or_case)

Get a `OrderedDict` containing all wells in the model or simulation case.
"""
function get_model_wells(model::MultiModel)
    wells = OrderedDict{Symbol, Any}()
    for (k, m) in pairs(model.models)
        if model_or_domain_is_well(m)
            wells[k] = physical_representation(m.data_domain)
        end
    end
    return wells
end

"""
    reservoir_system(flow = flow_system, thermal = thermal_system)

Set up a [`Jutul.CompositeSystem`](@ref) that combines multiple systems
together. In some terminologies this is refered to as a multi-physics system.
The classical example is to combine a flow system and a thermal system to create
a coupled flow and heat system.
"""
function reservoir_system(; flow = missing, thermal = missing, kwarg...)
    carg = Pair{Symbol, Jutul.JutulSystem}[]
    if !ismissing(flow)
        flow::MultiPhaseSystem
        push!(carg, :flow => flow)
    end
    if !ismissing(thermal)
        thermal::ThermalSystem
        push!(carg, :thermal => thermal)
    end
    for (k, v) in kwarg
        v::Jutul.JutulSystem
        push!(carg, k => v)
    end
    return CompositeSystem(:Reservoir; carg...)
end

"""
    model, parameters = setup_reservoir_model(reservoir, system; wells = [], <keyword arguments>)
    model, parameters = setup_reservoir_model(reservoir, system; wells = [w1, w2], backend = :csr, <keyword arguments>)

Set up a reservoir `MultiModel` for a given reservoir `DataDomain` typically set
up from  [`reservoir_domain`](@ref) and an optional vector of wells that are
created using [`setup_vertical_well`](@ref) and  [`setup_well`](@ref).

The routine automatically sets up a facility and couples the wells with the
reservoir and that facility.

# Keyword arguments

- `wells=[]`: Vector of wells (e.g. from [`setup_well`](@ref)) that are to be
  used in the model. Each well must have a unique name.
- `extra_out=true`: Return both the model and the parameters instead of just the
  model.
- `thermal = false`: Add additional equations for conservation of energy and
  temperature as a primary variable.
- `kgrad=nothing`: Type of spatial discretization to use:
    - `:tpfa` or `nothing` gives standard two-point flux approximation (TPFA) with hard-coded
      two-point assembly
    - `:tpfa_test` gives TPFA with specialized finite-volume assembly. Should be
      similar in performance to `:tpfa`, but does not make use of threads.
    - `:avgmpfa` gives a consistent linear MPFA scheme that is more accurate for
      meshes with anisotropic perm or non-orthogonal cells than `:tpfa`.
    - `:ntpfa` gives a consistent nonlinear MPFA scheme (nonlinear version of
      `:avgmpfa` that preserves monotonicity)
- `extra_outputs=Symbol[]`: Extra output variables for reservoir model. Defaults
  to "typical" values seen in reservoir simulation. Valid values: Vector of
  symbols to be output, `true` for all variables and `false` for the minimal set
  required to restart simulations (typically only the primary variables and mass
  of each component)

## Advanced model setup
Advanced options govern internals of the simulator, like type of automatic
differentation, how equations are linearized and so on. These should not impact
simulation results beyond what is allowed for the model tolerances, but can
impact simulation speed.

- `split_wells=false`: Add a facility model for each well instead of one
  facility model that controls all wells. This must be set to `true` if you want
  to use MPI or nonlinear domain decomposition.
- `backend=:csc`: Backend to use. Can be `:csc` for serial compressed sparse
  column CSC matrix, `:csr` for parallel compressed sparse row matrix. `:csr` is
  a bit faster and is recommended when using MPI, HYPRE or multiple threads.
  `:csc` uses the default Julia format and is interoperable with other Julia
  libraries.
- `context=DefaultContext()`: Context used for entire model. Not recommended to
  set up manually, use `backend` instead.
- `assemble_wells_together=true`: Assemble wells in a single big matrix rather
  than many small matrices.
- `block_backend=true`: Use block sparse representation. This is needed by the
  iterative solvers and corresponding preconditioners. Setting this to `false`
  will result in a direct solver being used. In addition, equations will be
  assembled in an order similar to that of MRST (equation major instead of cell
  major).
- `general_ad=false`: Use more general form of AD. Will result in slower
  execution speed than if set to true, but can be useful when working with
  custom discretizations.

## Increment and variable options
These options govern the range of values and the maximum allowable change of
properties over a single Newton iteration. Changing values for maximum change
will not change the equations themselves, but the values will change the rate of
nonlinear solver convergence. Typically, smaller values are more conservative
and reduce numerical difficulties, but can significantly increase the number of
iterations and the reduce the length of the average time-step. Setting very
small values can make it infeasible to solve the problems in a reasonable time.

Note that relative values are usually given relative to the cell value. If your
expected output values are close to zero (e.g. for near-atmospheric pressure)
low values can lead to slow convergence.

- `dp_max_abs=nothing`: Maximum allowable pressure change in SI units (Pa)
- `dp_max_rel=0.2`: Maximum allowable relative pressure change (default is 20%)
- `dp_max_abs_well=convert_to_si(50, :bar)`: Maximum allowable pressure change
  for wells in SI units (Pa)
- `dp_max_rel_well=nothing`: Maximum allowable relative pressure change in well
- `ds_max=0.2`: Maximum change in saturations
- `dz_max=0.2`: Maximum change in composition (for compositional models only)
- `p_min=JutulDarcy.DEFAULT_MINIMUM_PRESSURE`: Minimum pressure in model (hard limit)
- `p_max=Inf`: Maximum pressure in model (hard limit)
- `dr_max=Inf`: Maximum change in Rs/Rv for blackoil models over a Newton
  iteration. Taken relative to the saturated value of the cell.
- `dT_max_rel=nothing`: Maximum relative change in temperature (JutulDarcy uses Kelvin,
  so comments about changing limits near zero above does not apply to typical
  reservoir temperatures)
- `dT_max_abs=50.0`: Maximum absolute change in temperature (in °K/°C)
- `fast_flash=false`: Shorthand to enable `flash_reuse_guess` and
  `flash_stability_bypass`. These options can together speed up the time spent
  in flash solver for compositional models. Options are based on "Increasing the
  Computational Speed of Flash Calculations With Applications for Compositional,
  Transient Simulations" by Rasmussen et al (2006).
- `flash_reuse_guess=fast_flash`: Reuse previous flash guess when a cell remains
  in two-phase.
- `flash_stability_bypass=fast_flash`: Bypass stability testing for cells
  outside the two-phase and shadow region.
- `can_shut_wells=true`: Configure facility to allow shutting wells that
  repeatedly get rates with the wrong side. Disabling this can make certain
  models infeasible to simulate, but it can be useful to do so for simple models
  where you know that the wells should be operational.

"""
function setup_reservoir_model(reservoir::DataDomain, system::JutulSystem;
        wells = [],
        context = DefaultContext(),
        reservoir_context = nothing,
        general_ad = false,
        backend = :csc,
        thermal = false,
        extra_outputs = [:LiquidMassFractions, :VaporMassFractions, :Rs, :Rv, :Saturations],
        split_wells = false,
        assemble_wells_together = true,
        extra_out = true,
        dp_max_abs = nothing,
        dp_max_rel = 0.2,
        dp_max_abs_well = convert_to_si(50, :bar),
        dp_max_rel_well = nothing,
        ds_max = 0.2,
        dz_max = 0.2,
        p_min = DEFAULT_MINIMUM_PRESSURE,
        p_max = Inf,
        dr_max = Inf,
        dT_max_rel = nothing,
        dT_max_abs = 50.0,
        fast_flash = false,
        can_shut_wells = true,
        flash_reuse_guess = fast_flash,
        flash_stability_bypass = fast_flash,
        parameters = Dict{Symbol, Any}(),
        block_backend = true,
        nthreads = Threads.nthreads(),
        minbatch = 1000,
        kgrad = nothing,
        immutable_model = false,
        wells_systems = missing
    )
    if !(wells isa AbstractArray)
        wells = [wells]
    end
    # List of models (order matters)
    models = OrderedDict{Symbol, Jutul.AbstractSimulationModel}()
    reservoir_context, context = Jutul.select_contexts(
        backend; 
        main_context = reservoir_context,
        context = context,
        block_backend = block_backend,
        nthreads = nthreads,
        minbatch = minbatch
    )
    # We first set up the reservoir
    rmodel = SimulationModel(
        reservoir,
        system;
        context = reservoir_context,
        general_ad = general_ad,
        kgrad = kgrad
    )
    if thermal
        rmodel = add_thermal_to_model!(rmodel)
    end
    set_discretization_variables!(rmodel)
    set_reservoir_variable_defaults!(rmodel,
        dp_max_abs = dp_max_abs,
        dp_max_rel = dp_max_rel,
        p_min = p_min,
        p_max = p_max,
        dr_max = dr_max,
        ds_max = ds_max,
        dz_max = dz_max,
        dT_max_rel = dT_max_rel,
        dT_max_abs = dT_max_abs,
        flash_reuse_guess = flash_reuse_guess,
        flash_stability_bypass = flash_stability_bypass
    )

    if extra_outputs isa Bool
        if extra_outputs
            for k in keys(rmodel.secondary_variables)
                push!(rmodel.output_variables, k)
            end
            unique!(rmodel.output_variables)
        end
    else
        for k in extra_outputs
            if haskey(rmodel.secondary_variables, k)
                push!(rmodel.output_variables, k)
            end
            unique!(rmodel.output_variables)
        end
    end
    if haskey(reservoir, :diffusion) || haskey(reservoir, :diffusivity)
        rmodel.parameters[:Diffusivities] = Diffusivities()
    end
    models[:Reservoir] = rmodel
    # Then we set up all the wells
    mode = PredictionMode()
    if length(wells) > 0
        for (well_no, w) in enumerate(wells)
            if ismissing(wells_systems)
                wsys = system
            else
                wsys = wells_systems[well_no]
            end
            w_domain = DataDomain(w)
            wc = w.perforations.reservoir
            c = map_well_nodes_to_reservoir_cells(w, reservoir)
            for propk in [:temperature, :pvtnum]
                if haskey(reservoir, propk)
                    w_domain[propk] = reservoir[propk][c]
                end
            end
            wname = w.name
            wmodel = SimulationModel(w_domain, system, context = context)
            if thermal
                wmodel = add_thermal_to_model!(wmodel)
            end
            set_reservoir_variable_defaults!(wmodel,
                dp_max_abs = dp_max_abs_well,
                dp_max_rel = dp_max_rel_well,
                p_min = p_min,
                p_max = p_max,
                dr_max = dr_max,
                ds_max = ds_max,
                dz_max = dz_max,
                dT_max_rel = dT_max_rel,
                dT_max_abs = dT_max_abs,
                flash_reuse_guess = flash_reuse_guess,
                flash_stability_bypass = flash_stability_bypass
            )
            models[wname] = wmodel
            if split_wells
                wg = WellGroup([wname], can_shut_wells = can_shut_wells)
                F = SimulationModel(wg, mode, context = context, data_domain = DataDomain(wg))
                models[Symbol(string(wname)*string(:_ctrl))] = F
            end
        end
        # Add facility that groups the wells
        if !split_wells
            wg = WellGroup(map(x -> x.name, wells), can_shut_wells = can_shut_wells)
            F = SimulationModel(wg, mode, context = context, data_domain = DataDomain(wg))
            models[:Facility] = F
        end
    end

    # Put it all together as multimodel
    model = reservoir_multimodel(models, split_wells = split_wells, assemble_wells_together = assemble_wells_together, immutable_model = immutable_model)
    if extra_out
        parameters = setup_parameters(model, parameters)
        out = (model, parameters)
    else
        out = model
    end
    return out
end

function setup_reservoir_model(reservoir::DataDomain, label::Symbol; kwarg...)
    return setup_reservoir_model(reservoir, Val(label); kwarg...)
end

function set_reservoir_variable_defaults!(model;
        p_min,
        p_max,
        dp_max_abs,
        dp_max_rel,
        ds_max,
        dz_max,
        dr_max,
        dT_max_rel = nothing,
        dT_max_abs = nothing,
        flash_reuse_guess = false,
        flash_stability_bypass = flash_reuse_guess
    )
    # Replace various variables - if they are available
    replace_variables!(model, OverallMoleFractions = OverallMoleFractions(dz_max = dz_max), throw = false)
    replace_variables!(model, Saturations = Saturations(ds_max = ds_max), throw = false)
    replace_variables!(model, Temperature = Temperature(max_rel = dT_max_rel, max_abs = dT_max_abs), throw = false)
    replace_variables!(model, ImmiscibleSaturation = ImmiscibleSaturation(ds_max = ds_max), throw = false)
    replace_variables!(model, BlackOilUnknown = BlackOilUnknown(ds_max = ds_max, dr_max = dr_max), throw = false)

    rmodel = reservoir_model(model)
    if rmodel isa CompositionalModel
        sys = rmodel.system
        flash = FlashResults(rmodel, stability_bypass = flash_stability_bypass, reuse_guess = flash_reuse_guess)
        replace_variables!(rmodel, FlashResults = flash, throw = false)
    end
    p_def = Pressure(
        max_abs = dp_max_abs,
        max_rel = dp_max_rel,
        minimum = p_min,
        maximum = p_max
    )
    replace_variables!(model, Pressure = p_def, throw = false)
    return model
end

"""
    setup_reservoir_simulator(models, initializer, parameters = nothing; <keyword arguments>)

# Arguments
- `models`: either a single model or a Dict with the key :Reservoir for multimodels
- `initializer`: used to setup state0, must be compatible with `model`
- `parameters`: initialized parameters, must be compatible with `model` if provided

# Keyword arguments
- `split_wells`: Add facility model to each well (needed for domain decomposition and MPI solves)
- `assemble_wells_together`: Option to split wells into multiple sparse matrices (false argument experimental)
- `specialize=false`: use deep specialization of storage for faster execution, but significantly more compile time

Additional keyword arguments are documented in the version of this function that uses `JutulCase` as the input.
"""
function setup_reservoir_simulator(models, initializer, parameters = nothing;
        specialize = false,
        split_wells = false,
        assemble_wells_together = true,
        kwarg...
    )
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

"""
    setup_reservoir_simulator(case::JutulCase; <keyword arguments>)

# Keyword arguments

- `mode=:default`: Mode used for solving. Can be set to `:mpi` if running in MPI
  mode together with HYPRE, PartitionedArrays and MPI in your environment.
- `method=:newton`: Can be `:newton`, `:nldd` or `:aspen`. Newton is the most
  tested approach and `:nldd` can speed up difficult models. The `:nldd` option
  enables a host of additional options (look at the simulator config for more
  details).
- `presolve_wells=false`: Solve wells before each linearization. Can improve
  convergence for models with multisegment wells.

## Linear solver options

- `linear_solver=:bicgstab`: iterative solver to use (provided model supports
  it). Typical options are `:bicgstab` or `:gmres` Can alternatively pass a
  linear solver instance.
- `precond=:cpr`: preconditioner for iterative solver: Either :cpr or :ilu0.
- `rtol=1e-3`: relative tolerance for linear solver
- `linear_solver_arg`: `Dict` containing additional linear solver arguments.

## Timestepping options

- `initial_dt=si_unit(:day)`: initial timestep in seconds (one day by default)
- `target_ds=Inf`: target saturation change over a timestep used by timestepper.
- `target_its=8`: target number of nonlinear iterations per time step
- `offset_its=1`: dampening parameter for time step selector where larger values
  lead to more pessimistic estimates.
- `timesteps=:auto`: Set to `:auto` to use automatic timestepping, `:none` for
  no automatic timestepping (i.e. try to solve exact report steps)
- `max_timestep=si_unit(:year)`: Maximum internal timestep used in solver.

## Convergence criterions
- `tol_cnv=1e-3`: maximum allowable point-wise error (volume-balance)
- `tol_mb=1e-7`: maximum alllowable integrated error (mass-balance)
- `tol_cnv_well=10*tol_cnv`: maximum allowable point-wise error for well node
  (volume-balance)
- `tol_mb_well=1e4*tol_mb`: maximum alllowable integrated error for well node
  (mass-balance)
- `inc_tol_dp_abs=Inf`: Maximum allowable pressure change (absolute)
- `inc_tol_dp_rel=Inf`: Maximum allowable pressure change (absolute)
- `inc_tol_dz=Inf`: Maximum allowable composition change (compositional only).

## Inherited keyword arguments

Additional keyword arguments come from the base Jutul simulation framework. We
list a few of the most relevant entries here for convenience:
- `info_level = 0`: Output level. Set to 0 for minimal output, -1 for no output
  and 1 or more for increasing verbosity.
- `output_path`: Path to write output to.
- `max_nonlinear_iterations=15`: Maximum Newton iterations before a time-step is
  cut.
- `min_nonlinear_iterations=1`: Minimum number of Newtons to perform before
  checking convergence.
- `relaxation=Jutul.NoRelaxation()`: Dampening used for solves. Can be set to
  `Jutul.SimpleRelaxation()` for difficult models. Equivialent option is to set
  `true` for relaxation and `false` for no relaxation.
- `failure_cuts_timestep=true`: Cut timestep instead of throwing an error when
  numerical issues are encountered (e.g. linear solver divergence).
- `max_timestep_cuts=25`: Maximum number of timestep cuts before a solver gives
  up. Note that when using dynamic timestepping, this in practice defines a
  minimal timestep, with more than the prescribed number of cuts being allowed
  if the timestep is dynamically increased after cutting.
- `timestep_max_increase=10.0`: Max allowable factor to increase time-step by.
  Overrides any choices made in dynamic step selection.
- `timestep_max_decrease=0.1`: Max allowable factor to decrease time-step by.
  Overrides any choices made in dynamic step selection.
- `tol_factor_final_iteration=1.0`: If set to a value larger than 1.0, the final
  convergence check before a time-step is cut is relaxed by multiplying all
  tolerances with this value. Warning: Setting it to a large value can have
  severe impact on numerical accuracy. A value of 1 to 10 is typically safe if
  your default tolerances are strict.

"""
function setup_reservoir_simulator(case::JutulCase;
        mode = :default,
        method = :newton,
        precond = :cpr,
        linear_solver = :bicgstab,
        max_timestep = si_unit(:year),
        max_dt = max_timestep,
        rtol = nothing,
        initial_dt = si_unit(:day),
        target_ds = Inf,
        target_its = 8,
        offset_its = 1,
        tol_cnv = 1e-3,
        tol_mb = 1e-7,
        info_level = 0,
        tol_cnv_well = 10*tol_cnv,
        tol_mb_well = 1e4*tol_mb,
        tol_dp_well = 1e-3,
        inc_tol_dp_abs = Inf,
        inc_tol_dp_rel = Inf,
        failure_cuts_timestep = true,
        max_timestep_cuts = 25,
        inc_tol_dz = Inf,
        set_linear_solver = true,
        timesteps = :auto,
        relaxation = false,
        presolve_wells = false,
        parray_arg = Dict{Symbol, Any}(),
        linear_solver_arg = Dict{Symbol, Any}(),
        extra_timing_setup = false,
        nldd_partition = missing,
        nldd_arg = Dict{Symbol, Any}(),
        kwarg...
    )
    set_linear_solver = set_linear_solver || linear_solver isa Symbol
    # Handle old kwarg...
    max_timestep = min(max_dt, max_timestep)
    extra_kwarg = Dict{Symbol, Any}()
    # Setup simulator
    sim_kwarg = Dict{Symbol, Any}()
    sim_kwarg[:extra_timing] = extra_timing_setup
    if presolve_wells
        sim_kwarg[:prepare_step_handler] = PrepareStepWellSolver()
    end
    if mode == :default
        # Single-process solve
        if method == :newton
            sim = Simulator(case; sim_kwarg...)
        else
            extra_kwarg[:method] = method
            sim = NLDD.NLDDSimulator(case, nldd_partition; nldd_arg..., sim_kwarg...)
        end
    else
        # MPI/PArray solve
        if method == :newton
            make_sim = (m; kwarg...) -> Simulator(m; sim_kwarg..., kwarg...)
            pbuffer = false
        else
            make_sim = (m; kwarg...) -> NLDD.NLDDSimulator(m; nldd_arg..., sim_kwarg..., kwarg...)
            pbuffer = true
        end
        b = mode_to_backend(mode)
        sim = setup_reservoir_simulator_parray(case, b;
            simulator_constructor = make_sim,
            primary_buffer = pbuffer,
            parray_arg...
        )
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
        @assert isnothing(timesteps) || timesteps == :none
    end
    # Config: Linear solver, timestep selection defaults, etc...
    if set_linear_solver
        if info_level < 1
            v = -1
        else
            v = 0
        end
        if linear_solver isa Symbol
            extra_ls = (solver = linear_solver,)
        else
            extra_ls = NamedTuple()
        end
        extra_kwarg[:linear_solver] = reservoir_linsolve(case.model, precond;
            rtol = rtol,
            extra_ls...,
            linear_solver_arg...,
        )
    elseif isnothing(linear_solver)
        # Nothing
    else
        extra_kwarg[:linear_solver] = linear_solver
    end
    if relaxation isa Bool
        if relaxation
            relaxation = SimpleRelaxation()
        else
            relaxation = NoRelaxation()
        end
    end
    cfg = simulator_config(sim;
        extra_kwarg...,
        timestep_selectors = sel,
        info_level = info_level,
        relaxation = relaxation,
        max_timestep = max_timestep,
        failure_cuts_timestep = failure_cuts_timestep,
        max_timestep_cuts = max_timestep_cuts,
        kwarg...
    )
    set_default_cnv_mb!(cfg, sim,
        tol_cnv = tol_cnv,
        tol_mb = tol_mb,
        tol_cnv_well = tol_cnv_well,
        tol_mb_well = tol_mb_well,
        tol_dp_well = tol_dp_well,
        inc_tol_dp_abs = inc_tol_dp_abs,
        inc_tol_dp_rel = inc_tol_dp_rel,
        inc_tol_dz = inc_tol_dz
        )
    return (sim, cfg)
end


"""
    simulate_reservoir(state0, model, dt;
        parameters = setup_parameters(model),
        restart = false,
        forces = setup_forces(model),
        kwarg...
    )
    simulate_reservoir(case;
        kwarg...
    )

Convenience function for simulating a reservoir model. This function internally
calls [`setup_reservoir_simulator`](@ref), simulates the problem and returns a
[`ReservoirSimResult`](@ref). Keyword arguments are passed onto
[`setup_reservoir_simulator`](@ref) and are documented in that function.

You can optionally unpack this result into the most typical desired outputs:

`wellsols, states = simulate_reservoir(...)`

where `wellsols` contains the well results and `states` the reservoir results
(pressure, saturations and so on, in each cell of the reservoir domain).

# Examples
You can restart/resume simulations by both providing the `output_path` argument
and the `restart` argument:
```
# Automatically restart from last solved step and returning the outputs if the simulation was already solved.
result = simulate_reservoir(state0, model, dt, output_path = "/some/path", restart = true)

# Restart from step 5
result = simulate_reservoir(state0, model, dt, output_path = "/some/path", restart = 5)

# Start from the beginning (default)
result = simulate_reservoir(state0, model, dt, output_path = "/some/path", restart = false)
```

"""
function simulate_reservoir(state0, model, dt;
        parameters = setup_parameters(model),
        forces = setup_forces(model),
        kwarg...
    )
    case = JutulCase(model, dt, forces, state0 = state0, parameters = parameters)
    return simulate_reservoir(case; kwarg...)
end

function simulate_reservoir(case::JutulCase;
        config = missing,
        restart = false,
        simulator = missing,
        kwarg...
    )
    (; model, forces, state0, parameters, dt) = case
    if ismissing(simulator)
        sim, config_new = setup_reservoir_simulator(model, state0, parameters; kwarg...)
        if ismissing(config)
            config = config_new
        end
    else
        sim = simulator
        # May have been passed kwarg that should be accounted for
        if length(kwarg) > 0
            config = copy(config)
            for (k, v) in kwarg
                config[k] = v
            end
        end
        @assert !ismissing(config) "If simulator is provided, config must also be provided"
    end
    result = simulate!(sim, dt, forces = forces, config = config, restart = restart);
    return ReservoirSimResult(model, result, forces; simulator = sim, config = config)
end

function set_default_cnv_mb!(cfg::JutulConfig, sim::JutulSimulator; kwarg...)
    set_default_cnv_mb!(cfg, sim.model; kwarg...)
end

function set_default_cnv_mb!(cfg, model; kwarg...)
    set_default_cnv_mb_inner!(cfg[:tolerances], model; kwarg...)
end

function set_default_cnv_mb_inner!(tol, model;
        tol_cnv = 1e-3,
        tol_mb = 1e-7,
        tol_mb_well = 1e-3,
        tol_cnv_well = 1e-2,
        tol_dp_well = 1e-3,
        inc_tol_dp_abs = Inf,
        inc_tol_dp_rel = Inf,
        inc_tol_dz = Inf
        )
    sys = model.system
    if model isa Jutul.CompositeModel && hasproperty(model.system.systems, :flow)
        sys = flow_system(model.system)
    end
    if sys isa ImmiscibleSystem || sys isa BlackOilSystem || sys isa CompositionalSystem
        is_well = model_or_domain_is_well(model)
        if is_well
            if physical_representation(model) isa SimpleWell
                m = Inf
            else
                tol[:potential_balance] = (AbsMax = tol_dp_well,)
                m = tol_mb_well
            end
            c = tol_cnv_well
        else
            c = tol_cnv
            m = tol_mb
        end
        tol[:mass_conservation] = (
            CNV = c,
            MB = m,
            increment_dp_abs = inc_tol_dp_abs,
            increment_dp_rel = inc_tol_dp_rel,
            increment_dz = inc_tol_dz
        )
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

    has_flow = rmodel.system isa MultiPhaseSystem
    has_thermal = haskey(rmodel.equations, :energy_conservation)
    conservation = :mass_conservation
    energy = :energy_conservation
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

function reservoir_multimodel(models::AbstractDict;
        specialize = false,
        split_wells = false,
        immutable_model = false,
        assemble_wells_together = haskey(models, :Facility)
    )
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
    if immutable_model
        model = convert_to_immutable_storage(model)
    end
    return model
end

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
                init_w[:TotalMassFlux] = 0.0
            end
            c = map_well_nodes_to_reservoir_cells(wg, rmodel.data_domain)
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
        if eltype(v)<:AbstractFloat
            if !all(isfinite, v)
                jutul_message("setup_reservoir_state", "Non-finite entries found in initializer for $k.", color = :red)
            end
        end
        if isnothing(I)
            if !(k in svars)
                jutul_message("setup_reservoir_state", "Recieved primary variable $k, but this is not known to reservoir model.")
            end
        else
            push!(found, k)
        end
        res_init[k] = v
    end
    handle_alternate_primary_variable_spec!(res_init, found, rmodel, rmodel.system)
    if length(found) != length(pvars)
        missing_primary_variables = setdiff(pvars, found)
        @warn "Not all primary variables were initialized for reservoir model." missing_primary_variables
    end
    return setup_state(rmodel, res_init)
end

function handle_alternate_primary_variable_spec!(res_init, found, rmodel, system)
    # Internal utility to handle non-trivial specification of primary variables
    return res_init
end

"""
    setup_reservoir_forces(model; control = nothing, limits = nothing, set_default_limits = true, <keyword arguments>)

Set up driving forces for a reservoir model with wells
"""
function setup_reservoir_forces(model::MultiModel;
        control = nothing,
        limits = nothing,
        set_default_limits = true,
        bc = nothing,
        sources = nothing,
        kwarg...
    )
    submodels = model.models
    has_facility = any(x -> isa(x.domain, WellGroup), values(submodels))
    no_well_controls = isnothing(control) && isnothing(limits)
    if !isnothing(bc) && !(bc isa AbstractVector)
        bc = [bc]
    end
    if !isnothing(sources) && !(sources isa AbstractVector)
        sources = [sources]
    end
    reservoir_forces = (bc = bc, sources = sources)
    @assert no_well_controls || has_facility "Model must have facility when well controls are provided."
    if haskey(submodels, :Facility)
        # Unified facility for all wells
        facility = model.models.Facility

        surface_forces = setup_forces(facility,
            control = control,
            limits = limits,
            set_default_limits = set_default_limits
        )
        # Set up forces for the whole model.
        out = setup_forces(model, Facility = surface_forces; kwarg..., Reservoir = reservoir_forces)
    else
        new_forces = Dict{Symbol, Any}()
        for (k, m) in pairs(submodels)
            if model_or_domain_is_well(m) && !isnothing(control)
                ctrl_symbol = Symbol("$(k)_ctrl")
                @assert haskey(submodels, ctrl_symbol) "Controller for well $k must be present with the name $ctrl_symbol"
                subctrl = Dict{Symbol, Any}()
                subctrl[k] = control[k]
                sublimits = Dict{Symbol, Any}()
                if !isnothing(limits) && haskey(limits, k)
                    sublimits[k] = limits[k]
                else
                    sublimits[k] = nothing
                end
                facility = submodels[ctrl_symbol]
                new_forces[ctrl_symbol] = setup_forces(facility,
                    control = subctrl,
                    limits = sublimits,
                    set_default_limits = set_default_limits
                )
            end
        end
        out = setup_forces(model; pairs(new_forces)..., kwarg..., Reservoir = reservoir_forces)
    end
    # If the model is a composite model we need to do some extra work to pass on
    # flow forces with the correct label.
    #
    # TODO: At the moment we have no mechanism for setting up forces for thermal
    # specifically.
    for (k, m) in pairs(submodels)
        f = out[k]
        if m isa Jutul.CompositeModel
            mkeys = keys(m.system.systems)
            if haskey(f, :flow) && haskey(f, :thermal)
                tmp = f
            else
                tmp = Dict{Symbol, Any}()
                for mk in mkeys
                    tmp[mk] = nothing
                end
            end
            out[k] = (; pairs(tmp)...)
        end
    end
    return out
end


"""
    full_well_outputs(model, states, forces; targets = available_well_targets(model.models.Reservoir))

Get the full set of well outputs after a simulation has occured, for plotting or other post-processing.
"""
function full_well_outputs(model, states, forces; targets = missing)
    rmodel = reservoir_model(model)
    if ismissing(targets)
        targets = available_well_targets(rmodel)
    end
    cnames = component_names(rmodel.system)
    has_temperature = length(states) > 0 && haskey(first(states)[:Reservoir], :Temperature)
    well_t = Dict{Symbol, Union{Vector{Float64}, Vector{Symbol}}}
    out = Dict{Symbol, well_t}()
    for w in well_symbols(model)
        outw = well_t()
        for t in targets
            outw[translate_target_to_symbol(t(1.0))] = well_output(model, states, w, forces, t)
        end
        outw[:mass_rate] = well_output(model, states, w, forces, :TotalSurfaceMassRate)
        outw[:control] = well_output(model, states, w, forces, :control)
        if has_temperature
            outw[:temperature] = map(s -> s[w][:Temperature][well_top_node()], states)
        end
        for (i, cname) in enumerate(cnames)
            outw[Symbol("$(cname)_mass_rate")] = well_output(model, states, w, forces, i)
        end
        if haskey(outw, :lrat) && haskey(outw, :wrat)
            outw[:wcut] = outw[:wrat]./outw[:lrat]
        end
        if haskey(outw, :orat) && haskey(outw, :grat)
            outw[:gor] = abs.(outw[:grat])./max.(abs.(outw[:orat]), 1e-12)
        end
        out[w] = outw
    end
    return out
end

"""
    well_output(model, states, well_symbol, forces, target = BottomHolePressureTarget)

Get a specific well output from a valid operational target once a simulation is completed an `states` are available.
"""
function well_output(model::MultiModel, states, well_symbol, forces, target = BottomHolePressureTarget)
    n = length(states)

    well_number = 1
    for (k, m) in pairs(model.models)
        if k == well_symbol
            break
        end
        if model_or_domain_is_well(m)
            well_number += 1
        end
    end
    well_model = model.models[well_symbol]
    rhoS_o = reference_densities(well_model.system)

    to_target(t::DataType) = t(1.0)
    to_target(t::Type) = t(1.0)
    to_target(t::Symbol) = t
    to_target(t::Int) = t

    target_limit = to_target(target)
    if target == :control
        d = Symbol[]
        for (i, state) = enumerate(states)
            ctrl_grp = Symbol("$(well_symbol)_ctrl")
            if haskey(state, ctrl_grp)
                cfg = state[ctrl_grp][:WellGroupConfiguration]
            else
                cfg = state[:Facility][:WellGroupConfiguration]
            end
            ctrl = cfg.operating_controls[well_symbol]
            if ctrl isa DisabledControl
                t = :disabled
            else
                t = translate_target_to_symbol(ctrl.target)
            end
            push!(d, t)
        end
    else
        d = zeros(n)
        for (i, state) = enumerate(states)
            well_state = state[well_symbol]
            well_state = convert_to_immutable_storage(well_state)
            ctrl_grp = Symbol("$(well_symbol)_ctrl")
            if haskey(state, ctrl_grp)
                q_t = only(state[ctrl_grp][:TotalSurfaceMassRate])
            else
                q_t = state[:Facility][:TotalSurfaceMassRate][well_number]
            end
            if forces isa AbstractVector
                force = forces[i]
            else
                force = forces
            end
            if haskey(force, :outer)
                force = force.outer
            end
            if haskey(force, :Facility)
                gforce = force[:Facility]
            else
                gforce = force[Symbol("$(well_symbol)_ctrl")]
            end
            control = gforce.control[well_symbol]
            if target == :TotalSurfaceMassRate
                d[i] = q_t
            elseif target isa Int
                # Shorthand for component mass rate
                if control isa InjectorControl
                    mix = control.injection_mixture[target]
                else
                    totmass = well_state[:TotalMasses][:, 1]
                    mix = totmass[target]/sum(totmass)
                end
                d[i] = q_t*mix
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
                    current_control = replace_target(control, BottomHolePressureTarget(1.0))
                    rhoS, S = surface_density_and_volume_fractions(well_state)
                    v = well_target_value(q_t, current_control, target_limit, well_model, well_state, rhoS, S)
                    d[i] = v
                end
            end
        end
    end
    return d
end

"""
    well_symbols(model::MultiModel)

Get the keys of a `MultiModel` models that correspond to well models.
"""
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
    phases = get_phases(flow_system(model.system))
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
    if conn == :unit
        T = ones(Int, length(trans))
    else
        trans = max.(trans, 1e-20)
        if conn == :trans
            T = copy(trans)
            T = length(T)*T./sum(T)
            T = Int.(ceil.(T))
        elseif conn == :logtrans
            T = log10.(trans)
            offset = maximum(abs, T) + 1
            @. T += offset
            T = Int.(ceil.(T))
        else
            error("conn must be one of :trans, :unit or :logtrans, was $conn")
        end
    end
    groups = Vector{Vector{Int}}()
    if model isa MultiModel
        for (k, m) in pairs(model.models)
            wg = physical_representation(m.domain)
            if wg isa WellDomain
                rcells = vec(Int.(wg.perforations.reservoir))
                push!(groups, rcells)
            end
        end
    end
    return (N, T, groups)
end

function reservoir_partition(model::MultiModel, p)
    !haskey(model.models, :Facility) || throw(ArgumentError("Cannot partition model if split_wells = false in setup_reservoir_model"))
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
        wells = sr.wells.wells
        wk = keys(wells)
        nw = length(wk)
        print(io, "\n  wells ($nw present):\n")
        if nw > 0
            for k in wk
                print(io, "    :$k\n")
            end
            print(io, "    Results per well:\n")
            print_keys("       ", wells[first(wk)])
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

"""
    reservoir_groups_for_printing(model::MultiModel)

Break down the standard reservoir simulator model into groups based on if they
are the reservoir, wells, facility or something else.
"""
function reservoir_groups_for_printing(model::MultiModel)
    groups = Vector{Tuple{Symbol, Vector{Symbol}}}()

    rgroup = (:Reservoir, [:Reservoir])
    wgroup = (:Wells, Symbol[])
    fgroup = (:Facility, Symbol[])

    push!(groups, rgroup)
    push!(groups, wgroup)
    push!(groups, fgroup)

    for (k, m) in pairs(model.models)
        d = m.domain
        if d isa Jutul.DiscretizedDomain
            d = m.domain.representation
        end
        if k == :Reservoir
            continue
        elseif d isa JutulDarcy.WellDomain
            push!(wgroup[2], k)
        elseif d isa JutulDarcy.WellControllerDomain
            push!(fgroup[2], k)
        else
            push!(groups, (k, [k]))
        end
    end
    return groups
end

"""
    reservoir_transmissibility(d::DataDomain)
    reservoir_transmissibility(d::DataDomain; version = :xyz)

Special transmissibility function for reservoir simulation that handles
additional complexity present in industry grade models such as fault
multipliers, net-to-gross etc. The input should be a `DataDomain` instance
returned from [`reservoir_domain`](@ref)

The keyword argument `version` can be `:xyz` for permeability tensors that are
aligned with coordinate directions or `:ijk` to interpreted the permeability as
a diagonal tensor aligned with the logical grid directions. The latter choice is
only meaningful for a diagonal tensor.
"""
function reservoir_transmissibility(d::DataDomain; version = :xyz)
    g = physical_representation(d)
    N = d[:neighbors]
    nc = number_of_cells(d)
    faces, facepos = get_facepos(N, nc)
    facesigns = Jutul.get_facesigns(N, faces, facepos, nc)

    if version == :ijk
        face_dir = get_ijk_face_dir(g, N)
    else
        face_dir = missing
    end
    T_hf = compute_half_face_trans(
        d[:cell_centroids],
        d[:face_centroids],
        d[:normals],
        d[:areas],
        d[:permeability],
        faces, facepos, facesigns,
        version = version,
        face_dir = face_dir
    )
    nf = number_of_faces(d)
    neg_count = 0
    for (i, T_hf_i) in enumerate(T_hf)
        neg_count += T_hf_i < 0
        T_hf[i] = abs(T_hf_i)
    end
    # We only warn for significant amounts of negative transmissibilities, since
    # a few negative values is normal for reservoir grids.
    if neg_count > 0.1*length(T_hf)
        tran_tot = length(T_hf)
        perc = round(100*neg_count/tran_tot, digits = 2)
        jutul_message("Transmissibility", "Replaced $neg_count negative half-transmissibilities (out of $tran_tot, $perc%) with their absolute value.")
    end
    bad_count = 0
    for (i, T_hf_i) in enumerate(T_hf)
        if T_hf_i isa AbstractFloat && !isfinite(T_hf_i)
            bad_count += 1
            T_hf[i] = 0.0
        end
    end
    if bad_count > 0
        tran_tot = length(T_hf)
        perc = round(100*bad_count/tran_tot, digits = 2)
        jutul_message("Transmissibility", "Replaced $bad_count non-finite half-transmissibilities (out of $tran_tot, $perc%) with zero.")
    end
    if haskey(d, :net_to_gross)
        # Net to gross applies to vertical trans only
        otag = get_mesh_entity_tag(g, Faces(), :orientation, throw = false)
        if !ismissing(otag)
            # Use tags if provided
            face_is_vertical = map(1:nf) do face
                return mesh_entity_has_tag(g, Faces(), :orientation, :vertical, face)
            end
        elseif g isa CartesianMesh
            # Cartesian mesh is simple
            k_index = map(c -> cell_ijk(g, c), 1:nc)
            face_is_vertical = map(1:nf) do face
                l, r = N[:, face]
                return k_index[l] == k_index[r]
            end
        else
            # Fall back to normals
            normals = d[:normals]
            face_is_vertical = map(1:nf) do face
                nx, ny, nz = normals[:, face]
                return abs(nz) < max(abs(nx), abs(ny))
            end
        end
        for (c, ntg) in enumerate(d[:net_to_gross])
            if ntg isa AbstractFloat && ntg ≈ 1.0
                continue
            end
            for fp in facepos[c]:(facepos[c+1]-1)
                if face_is_vertical[faces[fp]]
                    T_hf[fp] *= ntg
                end
            end
        end
    end
    T = compute_face_trans(T_hf, N)
    if haskey(d, :transmissibility_multiplier, Faces())
        tm = d[:transmissibility_multiplier]
        @. T *= tm
    end
    if haskey(d, :numerical_aquifers)
        aquifers = d[:numerical_aquifers]
        bnd_areas = d[:boundary_areas]
        bnd_centroids = d[:boundary_centroids]
        cell_centroids = d[:cell_centroids]
        D = size(cell_centroids, 1)
        point_t = SVector{D, eltype(cell_centroids)}

        if haskey(d, :net_to_gross)
            ntg = d[:net_to_gross]
        else
            ntg = ones(nc)
        end
        T, num_aquifer_faces = set_aquifer_transmissibilities!(
            T,
            g, d[:permeability], ntg,
            aquifers,
            reinterpret(point_t, cell_centroids),
            reinterpret(point_t, bnd_centroids),
            g.boundary_faces.neighbors,
            bnd_areas
        )
    else
        num_aquifer_faces = 0
    end

    if haskey(d, :nnc)
        nnc = d[:nnc]
        num_nnc = length(nnc)
        # Aquifers come at the end, and NNC are just before the aquifer faces.
        # TODO: Do this in a less brittle way, e.g. by tags.
        offset = nf - num_nnc - num_aquifer_faces
        for (i, ncon) in enumerate(nnc)
            T[i + offset] = ncon[7]
        end
    end
    return T
end

function set_aquifer_transmissibilities!(T, mesh, perm, ntg, aquifers, cell_centroids, bnd_centroids, bnd_neighbors, bnd_areas)
    num_aquifer_faces = 0
    # Connections to the reservoir
    for (aq_id, aquifer) in pairs(aquifers)
        aqprm = aquifer.aquifer_cells[1]
        aquifer_cell = aqprm.cell
        R = aqprm.length/2.0
        for (bface, face, opt, tmult) in zip(
                aquifer.boundary_faces,
                aquifer.added_faces,
                aquifer.trans_option,
                aquifer.boundary_transmult
            )
            area_reservoir = bnd_areas[bface]
            reservoir_cell = bnd_neighbors[bface]
            dist = norm(bnd_centroids[bface] - cell_centroids[reservoir_cell])
            num_aquifer_faces += 1
            is_vertical = mesh_entity_has_tag(mesh, BoundaryFaces(), :orientation, :vertical, bface)

            if mesh_entity_has_tag(mesh, BoundaryFaces(), :ijk_orientation, :j, bface)
                dir = 2
            elseif mesh_entity_has_tag(mesh, BoundaryFaces(), :ijk_orientation, :k, bface)
                dir = 3
            else
                dir = 1
            end
            if is_vertical
                ntg_face = ntg[reservoir_cell]
            else
                ntg_face = 1.0
            end
            T_reservoir = perm[dir, reservoir_cell]*area_reservoir*ntg_face/dist

            if opt == 0
                area_aquifer = aqprm.area
            else
                @assert opt == 1 "Option for aquifer transmissibility expected to be 1 or 0, was $opt"
                area_aquifer = area_reservoir
            end
            T_aquifer = area_aquifer*aqprm.permeability/R
            effective_trans = tmult/(1.0/T_reservoir + 1.0/T_aquifer)
            if isfinite(effective_trans)
                T[face] = effective_trans
            else
                @error "Non-finite aquifer transmissibility for numerical aquifer $aq_id, setting to zero" T_aquifer T_reservoir aqprm
                T[face] = 0.0
            end
        end
    end
    # Aquifer internal connections
    for (aq_id, aquifer) in pairs(aquifers)
        aqprms = aquifer.aquifer_cells
        @assert length(aquifer.aquifer_faces) == length(aqprms)-1
        for (i, face) in enumerate(aquifer.aquifer_faces)
            num_aquifer_faces += 1
            curr = aqprms[i]
            next = aqprms[i+1]
            T_c = curr.area*curr.permeability/(curr.length/2.0)
            T_n = next.area*next.permeability/(next.length/2.0)
            T[face] = 1.0/(1.0/T_c + 1.0/T_n)
        end
    end
    return (T, num_aquifer_faces)
end

function get_ijk_face_dir(g, N)
    nf = number_of_faces(g)
    face_dir = ones(Int, nf)
    ijk_c = map(i -> cell_ijk(g, i), 1:number_of_cells(g))
    gdim = dim(g)
    for i in 1:nf
        l, r = N[:, i]
        ijk_l = ijk_c[l]
        ijk_r = ijk_c[r]
        for j in 1:gdim
            if ijk_l[j] != ijk_r[j]
                face_dir[i] = j
                break
            end
        end
    end
    return face_dir
end

function reservoir_regions(m::MultiModel, arg...; kwarg...)
    rm = reservoir_model(m)
    return reservoir_regions(rm, arg...; kwarg...)
end

function reservoir_regions(m::SimulationModel, arg...; kwarg...)
    return reservoir_regions(m.data_domain, arg...; kwarg...)
end

"""
    reservoir_regions(d::DataDomain, type = :pvtnum)

Get regions stored in `DataDomain`, returning `nothing` if the field is not
found.
"""
function reservoir_regions(d::DataDomain, type = :pvtnum)
    if haskey(d, type, Cells())
        out = d[type]
    else
        out = nothing
    end
    return out
end

function set_discretization_variables!(model::MultiModel)
    set_discretization_variables!(reservoir_model(model))
    return model
end

function set_discretization_variables!(model; ntpfa_potential = true)
    disc = model.domain.discretizations
    flow = model.domain.discretizations.mass_flow
    if flow isa PotentialFlow
        if eltype(flow.kgrad) != TPFA
            has_pc = !isnothing(get_variable(model, :CapillaryPressure, throw = false))
            # Potential, z on faces
            xyz = model.data_domain[:cell_centroids]
            if size(xyz, 1) == 3
                z = xyz[3, :]
                zmin, zmax = extrema(z)
                has_gravity = abs(zmax - zmin) > 1e-10
            else
                has_gravity = false
            end
            if ntpfa_potential || has_gravity || has_pc
                pp = PhasePotentials()
                acd = AdjustedCellDepths()
                if model.system isa CompositeSystem
                    pp = Pair(:flow, pp)
                    acd = Pair(:flow, acd)
                end
                set_secondary_variables!(model, PhasePotentials = pp)
                set_parameters!(model, AdjustedCellDepths = acd)
            end
        end
    end
    return model
end

export co2_inventory

"""
    inventory = co2_inventory(model, ws, states, t; cells = missing, co2_name = "CO2")
    inventory = co2_inventory(model, result::ReservoirSimResult; cells = missing)

Compute CO2 inventory for each step for a given `model`, well results `ws` and
reporting times t. If provided, the keyword argument `cells` will compute
inventory inside the region defined by the cells, and let any additional CO2 be
categorized as "outside region".

The inventory will be a Vector of Dicts where each entry contains a breakdown of
the status of the CO2 at that time, including residual and dissolution trapping.
"""
function co2_inventory(model, ws, states, t; cells = missing, co2_name = "CO2")
    dt = diff(vcat(0, t))
    @assert all(dt .> 0.0) "Fourth argument t should be time from simulation start. Maybe you passed timesteps?"
    rmodel = reservoir_model(model)
    domain = reservoir_domain(rmodel)
    nc = number_of_cells(domain)
    if ismissing(cells)
        cells = 1:nc
    end
    # Map component to index and do some checking on system
    sys = rmodel.system
    # Figure out phases
    phases = get_phases(sys)
    @assert length(phases) == 2 "This routine is only implemented for two phases"
    vapor = findfirst(isequal(VaporPhase()), phases)
    @assert !isnothing(vapor) "Did not find VaporPhase in get_phases: $phases"
    aqua = findfirst(isequal(AqueousPhase()), phases)
    liq = findfirst(isequal(LiquidPhase()), phases)
    if isnothing(aqua)
        aqua = liq
    end
    @assert !isnothing(aqua)
    phasepair = (aqua, vapor)

    immiscible = sys isa ImmiscibleSystem
    if immiscible
        cnames = component_names(sys)
        co2_index = vapor
        co2_name = cnames[co2_index]
    else
        cnames = component_names(sys)
        co2_index = findfirst(isequal(co2_name), cnames)
        @assert !isnothing(co2_index) "Did not find $co2_name in component_names: "
    end
    rate_name = Symbol("$(co2_name)_mass_rate")
    # Accumulate up well results
    nstep = length(dt)
    injected = zeros(nstep)
    produced = zeros(nstep)
    for (wk, wres) in pairs(ws.wells)
        q_co2 = wres[rate_name]
        for (i, v) in enumerate(q_co2)
            if v > 0
                injected[i] += v
            else
                injected[i] += abs(v)
            end
        end
    end
    relperm = rmodel[:RelativePermeabilities]
    relperm::ReservoirRelativePermeabilities

    return map(
        i -> co2_inventory_for_step(states, dt, injected, produced, i, relperm, co2_name, co2_index, phasepair, cells),
        eachindex(states, dt)
    )
end

function co2_inventory(model, result::ReservoirSimResult; kwarg...)
    ws, states, t = result
    return co2_inventory(model, ws, states, t; kwarg...)
end

function co2_inventory_for_step(states, timesteps, injected, produced, step_no, relperm, co2_name, co2_c_index, phases, cells)
    liquid, vapor = phases
    is_hyst = hysteresis_is_active(relperm)
    state = states[step_no]
    dt = timesteps[step_no]
    # Sum up injected up to this point
    inj_tot = 0.0
    # Sum up produced up to this point
    prod_tot = 0.0
    for i in 1:step_no
        dt_i = timesteps[i]
        inj_tot += injected[i]*dt_i
        prod_tot += produced[i]*dt_i
    end

    S = state[:Saturations]
    rho = state[:PhaseMassDensities]
    immiscible = !haskey(state, :LiquidMassFractions) && !haskey(state, :VaporMassFractions)
    if immiscible
        X = Y = missing
    else
        X = state[:LiquidMassFractions]
        Y = state[:VaporMassFractions]
    end
    total_masses = state[:TotalMasses]

    krg = relperm.krg
    regions = relperm.regions

    residual_trapped, mobile, dissolution_trapped, mass_inside = sum_co2_inventory(total_masses, krg, regions, X, Y, S, rho, cells, is_hyst, immiscible, liquid, vapor, co2_c_index, step_no)

    total_co2_mass = sum(view(total_masses, co2_c_index, :))
    return Dict(
        :residual => residual_trapped,
        :mobile => mobile,
        :dissolved => dissolution_trapped,
        :inside_total => mass_inside,
        :outside_total => total_co2_mass - mass_inside,
        :domain_total => total_co2_mass,
        :outside_domain => inj_tot - total_co2_mass - prod_tot,
        :injected => inj_tot,
        :produced => prod_tot
    )
end

function sum_co2_inventory(total_masses, krg, regions, X, Y, S, rho, cells, is_hyst, immiscible, liquid, vapor, co2_c_index, step_no)
    mass_inside = 0.0
    mass_outside = 0.0
    residual_trapped = 0.0
    dissolution_trapped = 0.0
    mobile = 0.0

    bad_cells = Int[]
    for c in cells
        res, free, diss, total, val = co2_inventory_for_cell(total_masses, krg, regions, X, Y, S, rho, c, is_hyst, immiscible, liquid, vapor, co2_c_index)
        if total > 1e-3 && abs(val - total)/max(total, 1e-3) > 1e-3
            push!(bad_cells, c)
        end

        residual_trapped += res
        mobile += free
        dissolution_trapped += diss
        mass_inside += total
    end
    if length(bad_cells) > 0
        jutul_message("CO2 Inventory", "Inconsistent masses in $(length(bad_cells)) cells for step $step_no. Maybe tolerances were relaxed?", color = :yellow)
    end
    return (residual_trapped, mobile, dissolution_trapped, mass_inside)
end

function co2_inventory_for_cell(total_masses, krg, regions, X, Y, S, rho, c, is_hyst, immiscible, liquid, vapor, co2_c_index)
    total_mass = 0.0
    for i in axes(total_masses, 1)
        total_mass += total_masses[i, c]
    end
    reg = region(regions, c)
    if is_hyst
        krg_cell = imbibition_table_by_region(krg, reg)
    else
        krg_cell = table_by_region(krg, reg)
    end
    # Liquid values
    sl = S[liquid, c]
    rho_l = rho[liquid, c]

    # Gas values
    sg_crit = krg_cell.critical
    sg = S[vapor, c]
    rho_v = rho[vapor, c]
    if immiscible
        Y_co2 = 1.0
        X_co2 = 0.0
    else
        Y_co2 = Y[co2_c_index, c]
        X_co2 = X[co2_c_index, c]
    end

    # Estimate PV from total mass + 2 phase assumption
    pv = total_mass/(sg*rho_v + sl*rho_l)

    sg_r = min(sg, sg_crit)
    sg_m = sg - sg_r

    # Trapped but in vapor phase
    co2_density_in_vapor = rho_v*Y_co2*pv
    res = sg_r*co2_density_in_vapor
    free = sg_m*co2_density_in_vapor
    # Solubility
    diss = rho_l*X_co2*sl*pv
    # Total and check
    total = total_masses[co2_c_index, c]
    val = free + res + diss
    return (res, free, diss, total, val)
end

export generate_jutuldarcy_examples

"""
    generate_jutuldarcy_examples(
        pth = pwd(),
        name = "jutuldarcy_examples";
        makie = nothing,
        project = true,
        force = false
    )

Make a copy of all JutulDarcy examples in `pth` in a subfolder
`jutuldarcy_examples`. An error is thrown if the folder already exists, unless
the `force=true`, in which case the folder will be overwritten and the existing
contents will be permanently lost. If `project=true`, a `Project.toml` file will
be generated with the same dependencies as that of the doc build system that
contains everything needed to run the examples.

The `makie` argument allows replacing the calls to GLMakie in the examples with
another backend specified by a `String`. There are no checks performed if the
replacement is the name of a valid Makie backend.
"""
function generate_jutuldarcy_examples(
        pth = pwd(),
        name = "jutuldarcy_examples";
        makie = nothing,
        project = true,
        force = false
    )
    if !ispath(pth)
        error("Destionation $pth does not exist. Specify a folder.")
    end
    dest = joinpath(pth, name)
    jdir, = splitdir(pathof(JutulDarcy))
    ex_dir = joinpath(jdir, "..", "examples")

    if ispath(dest)
        if !force
            error("Folder $name already already exists in $pth. Specify force = true, or choose another name.")
        end
    end
    cp(ex_dir, dest, force = force)
    if !isnothing(makie)
        replace_makie_calls!(dest, makie)
    end
    proj_location = joinpath(jdir, "..", "docs", "Project.toml")
    cp(ex_dir, proj_location, force = true)
    chmod(ex_dir, 0o777, recursive = true)
    return dest
end

function replace_makie_calls!(dest, makie)
    makie::Union{Symbol, AbstractString}
    makie_str = "$makie"
    for (root, dirs, files) in walkdir(dest)
        for ex in files
            if !endswith(lowercase(ex), ".jl")
                # Don't mess with Project.toml
                continue
            end
            ex_pth = joinpath(root, ex)
            ex_lines = readlines(ex_pth, keep=true)
            f = open(ex_pth, "w")
            for line in ex_lines
                newline = replace(line, "GLMakie" => makie_str)
                print(f, newline)
            end
            close(f)
        end
        for dir in dirs
            replace_makie_calls!(joinpath(root, dir), makie)
        end
    end
end

export reservoir_measurables

function reservoir_measurables(case::JutulCase, ws, states; kwarg...)
    return reservoir_measurables(case.model, ws, states; kwarg...)
end

function reservoir_measurables(case::JutulCase, result::ReservoirSimResult; kwarg...)
    ws, states = result
    return reservoir_measurables(case.model, ws, states; kwarg...)
end

function reservoir_measurables(model, result::ReservoirSimResult; kwarg...)
    ws, states = result
    return reservoir_measurables(model, ws, states; kwarg...)
end

function reservoir_measurables(model, ws, states = missing; type = :field)
    @assert type == :field
    wells = ws.wells
    time = ws.time

    out = Dict{Symbol, Any}(:time => time)

    model = reservoir_model(model)
    reservoir = reservoir_domain(model)
    pv = pore_volume(reservoir)
    pv_t = sum(pv)
    sys = model.system
    is_blackoil = sys isa ImmiscibleSystem || sys isa StandardBlackOilSystem
    phases = get_phases(model.system)
    wix = findfirst(isequal(AqueousPhase()), phases)
    oix = findfirst(isequal(LiquidPhase()), phases)
    gix = findfirst(isequal(VaporPhase()), phases)

    has_water = !isnothing(wix)
    has_oil = !isnothing(oix)
    has_gas = !isnothing(gix)
    n = length(time)

    function add_entry(name, legend, unit = :id; is_rate = false)
        values = zeros(n)
        out[name] = (values = values, legend = legend, unit_type = unit, is_rate = is_rate)
        return values
    end
    # Production of different types
    flpr = add_entry(:flpr, "Field liquid production rate (oil + water)", :liquid_volume_surface, is_rate = true)
    fwpr = add_entry(:fwpr, "Field water production rate", :liquid_volume_surface, is_rate = true)
    fopr = add_entry(:fopr, "Field oil production rate", :liquid_volume_surface, is_rate = true)
    fgpr = add_entry(:fgpr, "Field gas production rate", :gas_volume_surface, is_rate = true)

    # Injection types
    fwir = add_entry(:fwir, "Field water injection rate", :liquid_volume_surface, is_rate = true)
    foir = add_entry(:foir, "Field oil injection rate", :liquid_volume_surface, is_rate = true)
    fgir = add_entry(:fgir, "Field gas injection rate", :gas_volume_surface, is_rate = true)

    # Reservoir values
    if is_blackoil
        fwip = add_entry(:fwip, "Field water component in place (surface volumes)", :liquid_volume_surface)
        foip = add_entry(:foip, "Field oil component in place (surface volumes)", :liquid_volume_surface)
        fgip = add_entry(:fgip, "Field gas component in place (surface volumes)", :gas_volume_surface)
    end

    fwipr = add_entry(:fwipr, "Field water in place (reservoir volumes)", :liquid_volume_reservoir)
    foipr = add_entry(:foipr, "Field oil in place (reservoir volumes)", :liquid_volume_reservoir)
    fgipr = add_entry(:fgipr, "Field gas in place (reservoir volumes)", :gas_volume_reservoir)

    fprh = add_entry(:fprh, "Field average pressure (hydrocarbon volume weighted)", :pressure)
    pres = add_entry(:pres, "Field average pressure", :pressure)

    function sum_well_rates!(vals, k::Symbol; is_prod::Bool)
        for (wk, wval) in pairs(wells)
            if is_prod
                for (i, v) in enumerate(wval[k])
                    vals[i] -= min(0.0, v)
                end
            else
                for (i, v) in enumerate(wval[k])
                    vals[i] += max(0.0, v)
                end
            end
        end
        return vals
    end
    if has_water
        sum_well_rates!(flpr, :wrat, is_prod = true)
        sum_well_rates!(fwpr, :wrat, is_prod = true)
        sum_well_rates!(fwir, :wrat, is_prod = false)
    end
    if has_oil
        sum_well_rates!(flpr, :orat, is_prod = true)
        sum_well_rates!(fopr, :orat, is_prod = true)
        sum_well_rates!(foir, :orat, is_prod = false)
    end
    if has_gas
        sum_well_rates!(fgpr, :grat, is_prod = true)
        sum_well_rates!(fgir, :grat, is_prod = false)
    end
    if !ismissing(states)
        if haskey(states[1], :Reservoir)
            states = map(x -> x[:Reservoir], states)
        end
        for (i, state) in enumerate(states)
            p = state[:Pressure]
            s = state[:Saturations]
            tm = state[:TotalMasses]
            if has_water
                fwipr[i] = sum(ix -> pv[ix]*s[wix, ix], eachindex(pv))
                if is_blackoil
                    fwip[i] = sum(tm[wix, :])
                end
            end
            if has_oil
                foipr[i] = sum(ix -> pv[ix]*s[oix, ix], eachindex(pv))
                if is_blackoil
                    foip[i] = sum(tm[oix, :])
                end
            end
            if has_gas
                fgipr[i] = sum(ix -> pv[ix]*s[gix, ix], eachindex(pv))
                if is_blackoil
                    fgip[i] = sum(tm[gix, :])
                end
            end

            mean_p = sum(ix -> pv[ix]*p[ix], eachindex(p, pv))/pv_t
            if has_water
                hc_p = sum(ix -> pv[ix]*p[ix]*(1.0 - s[wix, ix]), eachindex(p, pv))
                pv_hc = sum(ix -> pv[ix]*(1.0 - s[wix, ix]), eachindex(p, pv))
                hc_mean_p = hc_p/pv_hc
            else
                hc_mean_p = mean_p
            end
            fprh[i] = hc_mean_p
            pres[i] = mean_p
        end
    end
    # Derived quantities
    fwct = add_entry(:fwct, "Field production water cut")
    @. fwct = fwpr./flpr
    fgor = add_entry(:fgor, "Field gas-oil production ratio")
    @. fgor = fgpr./max.(fopr, 1e-12)

    return out
end
