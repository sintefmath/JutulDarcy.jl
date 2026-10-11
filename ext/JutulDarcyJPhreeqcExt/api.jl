function validate_solution_input(input)
    forbidden = ("EQUILIBRIUM_PHASES", "EXCHANGE", "SURFACE", "GAS_PHASE",
        "SOLID_SOLUTIONS", "SOLID_SOLUTION", "KINETICS", "REACTION", "MIX",
        "USE", "SAVE", "COPY", "DELETE", "RUN_CELLS", "TRANSPORT", "ADVECTION",
        "SELECTED_OUTPUT", "USER_PUNCH", "INCLUDE\$", "SOLUTION_RAW",
        "EQUILIBRIUM_PHASES_RAW", "REACTION_TEMPERATURE", "REACTION_PRESSURE")
    for line in split(input, '\n')
        text = strip(first(split(line, '#'; limit = 2)))
        isempty(text) && continue
        keyword = uppercase(first(split(text)))
        keyword in forbidden && throw(ArgumentError("$keyword is not supported in equilibrium solution input. Supply pure minerals through equilibrium_phases; include their inventory in tracer totals."))
    end
    solutions = UnitRange{Int}[]
    for match in eachmatch(r"(?im)^\s*SOLUTION\s+(\d+)(?:\s*-\s*(\d+))?", input)
        lo = parse(Int, match.captures[1])
        hi = isnothing(match.captures[2]) ? lo : parse(Int, match.captures[2])
        lo <= hi || throw(ArgumentError("Invalid SOLUTION range $lo-$hi."))
        push!(solutions, lo:hi)
    end
    any(r -> 1 in r, solutions) || throw(ArgumentError("Input must define SOLUTION 1 for initialization."))
    return solutions
end

function setup_phreeqc_model!(model::Union{SimulationModel, MultiModel}, input::AbstractString;
        equilibrium_phases = String[], database = "phreeqc.dat", temperature = 298.15,
        batch_size = 32, threads = Threads.nthreads(), difference = :forward,
        relative_step = 1e-5, concentration_scale = 1e-5,
        update_pore_volume = true, mineral_molar_volumes = nothing,
        reference_mineral_volume = 0.0, reuse_speciation = true)
    rmodel = reservoir_model(model)
    rmodel.context isa Jutul.KernelAbstractionsContext && throw(ArgumentError("PHREEQC chemistry requires a CPU context without secondary-variable kernels."))
    haskey(rmodel.equations, :tracers) && throw(ArgumentError("Add PHREEQC before setting up other tracers; existing tracer sets are not supported."))
    batch_size isa Integer && batch_size > 0 || throw(ArgumentError("batch_size must be a positive integer."))
    threads isa Integer && 0 < threads <= Threads.nthreads() || throw(ArgumentError("threads must be between 1 and Julia's thread count."))
    difference in (:forward, :central) || throw(ArgumentError("difference must be :forward or :central."))
    isfinite(relative_step) && relative_step > 0 || throw(ArgumentError("relative_step must be positive and finite."))
    isfinite(temperature) && temperature > 0 || throw(ArgumentError("temperature must be positive and in Kelvin."))
    isfinite(concentration_scale) && concentration_scale > 0 || throw(ArgumentError("concentration_scale must be positive and finite."))
    solutions = validate_solution_input(input)
    if model isa MultiModel
        for m in values(model.models)
            if JutulDarcy.model_or_domain_is_well(m)
                Jutul.physical_representation(m.domain) isa JutulDarcy.SimpleWell ||
                    throw(ArgumentError("PHREEQC currently supports simple wells only."))
                haskey(m.equations, :tracers) && throw(ArgumentError("A well already has tracers."))
            end
        end
    end
    phases = String[phase for phase in equilibrium_phases]
    length(unique(phases)) == length(phases) || throw(ArgumentError("Duplicate equilibrium phases."))
    all(p -> occursin(r"^[A-Za-z0-9_()+.:-]+$", p), phases) || throw(ArgumentError("Invalid equilibrium phase name."))
    update_pore_volume isa Bool || throw(ArgumentError("update_pore_volume must be a Bool."))
    reuse_speciation isa Bool || throw(ArgumentError("reuse_speciation must be a Bool."))
    if !isnothing(mineral_molar_volumes)
        length(mineral_molar_volumes) == length(phases) || throw(ArgumentError("Supply one mineral molar volume per equilibrium phase."))
        all(v -> isfinite(v) && v > 0, mineral_molar_volumes) || throw(ArgumentError("Mineral molar volumes must be positive, in m³/mol."))
    end
    phase_index = findfirst(p -> p isa AqueousPhase, JutulDarcy.get_phases(rmodel.system))
    isnothing(phase_index) && throw(ArgumentError("PHREEQC requires an aqueous flow phase."))
    nc = number_of_cells(rmodel.domain)
    nc > 0 || throw(ArgumentError("PHREEQC requires at least one cell."))
    reference_volume = reference_mineral_volume isa Real ? fill(Float64(reference_mineral_volume), nc) : Float64.(reference_mineral_volume)
    length(reference_volume) == nc && all(v -> isfinite(v) && v >= 0, reference_volume) ||
        throw(ArgumentError("reference_mineral_volume must contain a nonnegative volume in m³ per cell."))
    update_pore_volume && !isempty(phases) &&
        !any(vars -> haskey(vars, :FluidVolume), (rmodel.parameters, rmodel.secondary_variables)) &&
        throw(ArgumentError("Reactive pore volume requires a parameter or secondary FluidVolume."))
    batch_size = min(batch_size, nc)
    nworkers = min(threads, cld(nc, batch_size))
    db = isfile(database) ? abspath(database) : database_path(database)
    workers = ChemistryWorker[]
    components = String[]
    try
        for i in 1:nworkers
            worker, names = create_worker(String(input), phases, db, batch_size, temperature,
                mineral_molar_volumes, update_pore_volume, reuse_speciation)
            push!(workers, worker)
            if i == 1
                components = names
            else
                names == components || error("Inconsistent PHREEQC components between workers.")
            end
        end
        signed = [name in ("H", "O", "Charge") for name in components]
        reference = copy(first(workers).reference)
        molar_masses = zeros(length(components))
        native_read() do
            checked(RM_GetGfw(first(workers).id, molar_masses), first(workers).id, "component formula weights")
        end
        molar_masses ./= 1000.0
        pool = Channel{ChemistryWorker}(nworkers)
        foreach(w -> put!(pool, w), workers)
        n = length(components)
        np = length(phases)
        no = n + np + 1
        variable_temperature = any(vars -> haskey(vars, :Temperature),
            (rmodel.primary_variables, rmodel.secondary_variables, rmodel.parameters))
        lifecycle_lock = ReentrantLock()
        chemistry = PhreeqcEquilibrium(workers, pool, gensym(:PhreeqcWorker),
            components, molar_masses, phases, solutions, signed, reference, phase_index,
            Float64(temperature), variable_temperature, difference, Float64(relative_step),
            batch_size, zeros(n + 2, nc), zeros(no, nc), zeros(no, n + 2, nc),
            fill(false, nc), fill(false, nc), [ReentrantLock() for _ in 1:cld(nc, batch_size)],
            false, 0, lifecycle_lock, Threads.Condition(lifecycle_lock))
        tracers = [PhreeqcTracer(phase_index, Float64(concentration_scale)) for _ in components]
        add_tracers_to_model!(model, tracers, names = components)
        rmodel.secondary_variables[:PhreeqcPhaseSplit] = chemistry
        rmodel.secondary_variables[:PhreeqcAqueousConcentrations] = PhreeqcSplitConcentrations(1:n)
        if np > 0
            rmodel.secondary_variables[:PhreeqcMineralConcentrations] = PhreeqcSplitConcentrations((n + 1):(n + np))
            rmodel.secondary_variables[:PhreeqcMineralVolume] = PhreeqcMineralVolume(no)
            if update_pore_volume
                if haskey(rmodel.secondary_variables, :FluidVolume)
                    base_volume = :PhreeqcUnreactedFluidVolume
                    rmodel.secondary_variables[base_volume] = pop!(rmodel.secondary_variables, :FluidVolume)
                else
                    base_volume = :StaticFluidVolume
                    rmodel.parameters[base_volume] = pop!(rmodel.parameters, :FluidVolume)
                end
                rmodel.parameters[:PhreeqcReferenceMineralVolume] = PhreeqcReferenceMineralVolume(reference_volume)
                rmodel.secondary_variables[:FluidVolume] = PhreeqcPoreVolume(base_volume, phase_index)
            end
            append!(rmodel.output_variables, [:PhreeqcMineralConcentrations, :PhreeqcMineralVolume, :FluidVolume])
        end
        rmodel.secondary_variables[:PhreeqcImmobileConcentrations] = PhreeqcImmobileConcentrations(n)
        rmodel.primary_variables[:TracerConcentrations] = PhreeqcConcentrations(signed, reference)
        if model isa MultiModel
            for m in values(model.models)
                if JutulDarcy.model_or_domain_is_well(m)
                    m.primary_variables[:TracerConcentrations] = PhreeqcConcentrations(signed, reference)
                end
            end
        end
        append!(rmodel.output_variables, [:PhreeqcAqueousConcentrations, :PhreeqcImmobileConcentrations])
        unique!(rmodel.output_variables)
        finalizer(chemistry) do runtime
            # Finalizers must not block on locks held by an interrupted task.
            @async close_chemistry!(runtime)
        end
    catch
        for w in workers
            native_write() do
                RM_Destroy(w.id)
            end
        end
        rethrow()
    end
    return model
end

function chemistry_runtime(model)
    c = get(reservoir_model(model).secondary_variables, :PhreeqcPhaseSplit, nothing)
    c isa PhreeqcEquilibrium || throw(ArgumentError("Model has no PHREEQC chemistry."))
    return c
end

phreeqc_component_names(model::Union{SimulationModel, MultiModel}) = copy(chemistry_runtime(model).components)
phreeqc_component_molar_masses(model::Union{SimulationModel, MultiModel}) = copy(chemistry_runtime(model).molar_masses)
close_phreeqc!(model::Union{SimulationModel, MultiModel}) = close_chemistry!(chemistry_runtime(model))

function phreeqc_initial_concentrations(model::Union{SimulationModel, MultiModel}; solution = 1)
    c = chemistry_runtime(model)
    solution isa Integer && any(r -> solution in r, c.solutions) ||
        throw(ArgumentError("SOLUTION $solution is not defined in the chemistry input."))
    return lock(c.lifecycle_lock) do
        c.closed && error("The PHREEQC runtime has been closed.")
        vals = zeros(length(c.components))
        with_chemistry_worker(c) do worker
            native_read() do
                checked(RM_InitialPhreeqc2Concentrations(worker.id, vals, 1,
                    [Int(solution)], [-1], [1.0]), worker.id, "solution concentrations")
            end
        end
        vals .- c.reference
    end
end
