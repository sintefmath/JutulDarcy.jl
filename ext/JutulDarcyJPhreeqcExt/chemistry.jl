# PhreeqcRM's instance registry is shared by all modules. Serialize creation
# and destruction; workers only read the registry while evaluating chemistry.
const INSTANCE_LOCK = ReentrantLock()
const INSTANCE_CONDITION = Threads.Condition(INSTANCE_LOCK)
const ACTIVE_NATIVE_CALLS = Ref(0)

function native_read(f)
    lock(INSTANCE_LOCK) do
        ACTIVE_NATIVE_CALLS[] += 1
    end
    try
        return f()
    finally
        lock(INSTANCE_LOCK) do
            ACTIVE_NATIVE_CALLS[] -= 1
            notify(INSTANCE_CONDITION; all = true)
        end
    end
end

function native_write(f)
    return lock(INSTANCE_LOCK) do
        while ACTIVE_NATIVE_CALLS[] > 0
            wait(INSTANCE_CONDITION)
        end
        f()
    end
end

function checked(status, id, operation)
    if status < 0
        error("PHREEQC $operation failed ($status): $(RM_GetErrorString(id))")
    end
    return status
end

function create_worker(input, phases, database, batch_size, temperature, molar_volumes, update_pore_volume,
        reuse_speciation = true, native_threads = 1, rebalance_fraction = 0.0)
    return native_write() do
        create_worker_impl(input, phases, database, batch_size, temperature, molar_volumes,
            update_pore_volume, reuse_speciation, native_threads, rebalance_fraction)
    end
end

function create_worker_impl(input, phases, database, batch_size, temperature, molar_volumes,
        update_pore_volume, reuse_speciation, native_threads, rebalance_fraction)
    id = lock(INSTANCE_LOCK) do
        RM_Create(batch_size, native_threads)
    end
    id >= 0 || error("Could not create PhreeqcRM instance ($id).")
    try
        actual_threads = RM_GetThreadCount(id)
        actual_threads == native_threads || throw(ArgumentError(
            "PhreeqcRM provides $actual_threads threads, but $native_threads were requested. Native threading requires an OpenMP-enabled library."))
        checked(RM_SetErrorHandlerMode(id, 0), id, "error mode")
        checked(RM_SetScreenOn(id, 0), id, "screen output")
        checked(RM_SetErrorOn(id, 0), id, "error output")
        checked(RM_SetPrintChemistryOn(id, 0, 0, 0), id, "print output")
        checked(RM_SetSelectedOutputOn(id, 0), id, "selected output")
        checked(RM_SetSpeciesSaveOn(id, 0), id, "species storage")
        checked(RM_SetRebalanceFraction(id, rebalance_fraction), id, "native load balancing")
        checked(RM_SetComponentH2O(id, 0), id, "conserved H/O components")
        checked(RM_SetUnitsSolution(id, 3), id, "solution units")
        checked(RM_UseSolutionDensityVolume(id, 0), id, "transport volume")
        checked(RM_SetPartitionUZSolids(id, 0), id, "solid partitioning")
        checked(RM_SetPorosity(id, ones(batch_size)), id, "porosity")
        checked(RM_SetSaturationUser(id, ones(batch_size)), id, "saturation")
        checked(RM_SetRepresentativeVolume(id, ones(batch_size)), id, "volume")
        checked(RM_SetDensityUser(id, ones(batch_size)), id, "density")
        checked(RM_SetTemperature(id, fill(temperature - 273.15, batch_size)), id, "temperature")
        checked(RM_LoadDatabase(id, database), id, "database")
        checked(RM_RunString(id, 1, 1, 0,
            "KNOBS\n    -convergence_tolerance 1e-12\n    -iterations 200\n    -numerical_derivatives false\nEND\n"), id, "solver tolerances")
        # Thermodynamic definitions must reach workers as well as InitialPhreeqc.
        # InitialSolutions2Module below replaces their temporary solutions.
        # Mineral assemblages start at zero: all available mass is in totals.
        phase_text = isempty(phases) ? "" : "\nEQUILIBRIUM_PHASES 1\n" *
            join(("    $phase 0.0 0.0" for phase in phases), "\n") * "\nEND\n"
        checked(RM_RunString(id, 1, 1, 0, input * "\nEND\n"), id, "input")
        # PrintChemistryOn does not disable input echoing, status polling or
        # USER_GRAPH. Keep only the selected values required by the coupling.
        checked(RM_RunString(id, 1, 1, 0,
            "PRINT\n    -reset false\n    -echo_input false\n    -status false\n    -user_graph false\nEND\n"), id, "print controls")
        if !isempty(phase_text)
            checked(RM_RunString(id, 0, 1, 0, phase_text), id, "initial mineral definitions")
        end
        n = checked(RM_FindComponents(id), id, "components")
        components = [RM_GetComponent(id, i - 1) for i in 1:n]
        gfw = zeros(n)
        checked(RM_GetGfw(id, gfw), id, "component formula weights")
        h = findfirst(==("H"), components)
        o = findfirst(==("O"), components)
        reference = zeros(n)
        water_gfw = 2*gfw[h] + gfw[o]
        reference[h] = 2*gfw[h]/water_gfw
        reference[o] = gfw[o]/water_gfw
        checked(RM_InitialSolutions2Module(id, fill(1, batch_size)), id, "initial solutions")
        if !isempty(phases)
            checked(RM_InitialEquilibriumPhases2Module(id, fill(1, batch_size)), id, "initial minerals")
        end
        checked(RM_SetTimeStep(id, 0.0), id, "timestep")
        checked(RM_StateSave(id, 1), id, "baseline state")
        stoichiometry, volumes = mineral_properties(id, phases, components, gfw,
            batch_size, molar_volumes, update_pore_volume)
        np = length(phases)
        nd = n + 2
        no = n + np + 1
        worker = ChemistryWorker(id, zeros(nd, batch_size), zeros(nd, batch_size),
            zeros(n*batch_size), ones(batch_size), zeros(batch_size), reference,
            fill(1, batch_size), reuse_speciation, true,
            stoichiometry, volumes, zeros(2*np*batch_size), zeros(np, batch_size),
            zeros(no, batch_size), zeros(no, batch_size), zeros(no, nd, batch_size), zeros(batch_size))
        return worker, components
    catch
        lock(INSTANCE_LOCK) do
            RM_Destroy(id)
        end
        rethrow()
    end
end

# Only setup uses formula metadata. Runtime selected output contains the two
# quantities needed per mineral: amount and saturation index, with no files.
function mineral_properties(id, phases, components, gfw, nb, supplied_volumes, update_pore_volume)
    n, np = length(components), length(phases)
    stoichiometry = zeros(n, np)
    volumes = zeros(np)
    np == 0 && return stoichiometry, volumes
    lines = String[]
    number = 0
    function basic(statement)
        number += 10
        push!(lines, "$number $statement")
    end
    for phase in phases
        basic("formula\$ = PHASE_FORMULA\$(\"$phase\", count, elements\$, coefficients)")
        basic("PUNCH PHASE_VM(\"$phase\")")
        for component in components
            basic("amount = 0")
            basic("FOR k = 1 TO count")
            basic("IF elements\$(k) = \"$component\" THEN amount = coefficients(k)")
            basic("NEXT k")
            basic("PUNCH amount")
        end
    end
    headings = join(("metadata$i" for i in 1:np*(n + 1)), " ")
    definition = selected_output_definition(headings, join(lines, "\n"))
    checked(RM_SetSelectedOutputOn(id, 1), id, "mineral output")
    checked(RM_RunString(id, 1, 0, 0, definition), id, "mineral metadata")
    checked(RM_RunCells(id), id, "mineral metadata evaluation")
    RM_GetSelectedOutputColumnCount(id) == np*(n + 1) || error("Unexpected PHREEQC mineral metadata layout.")
    values = zeros(nb*np*(n + 1))
    checked(RM_GetSelectedOutput(id, values), id, "mineral metadata retrieval")
    for p in 1:np
        offset = (p - 1)*(n + 1)
        volumes[p] = isnothing(supplied_volumes) ? values[offset*nb + 1]*1e-6 : supplied_volumes[p]
        if update_pore_volume && !(isfinite(volumes[p]) && volumes[p] > 0)
            throw(ArgumentError("No positive molar volume for $(phases[p]). Supply mineral_molar_volumes in m³/mol or disable update_pore_volume."))
        end
        for i in 1:n
            stoichiometry[i, p] = values[(offset + i)*nb + 1]*gfw[i]/1000.0
        end
    end
    headings = join(("amount$p saturation$p" for p in 1:np), " ")
    punch = "10 PUNCH " * join(("EQUI(\"$phase\"), SI(\"$phase\")" for phase in phases), ", ")
    checked(RM_RunString(id, 1, 0, 0, selected_output_definition(headings, punch)), id, "mineral output definition")
    return stoichiometry, volumes
end

function selected_output_definition(headings, program)
    return "SELECTED_OUTPUT 1\n    -reset false\n    -high_precision true\n" *
        "USER_PUNCH 1\n    -headings $headings\n    -start\n$program\n    -end\nEND\n"
end

# A lease belongs to a task, and nested property calls find the same prepared
# instance. Scoped TLS avoids retaining a returned worker in long-lived tasks.
function with_chemistry_worker(f, chemistry)
    worker = get(task_local_storage(), chemistry.task_key, nothing)
    if !isnothing(worker)
        return f(worker)
    end
    worker = take!(chemistry.available)
    try
        return task_local_storage(() -> f(worker), chemistry.task_key, worker)
    finally
        put!(chemistry.available, worker)
    end
end

chemistry_worker(chemistry) = task_local_storage(chemistry.task_key)::ChemistryWorker

function equilibrium!(output, w::ChemistryWorker, x)
    isempty(w.molar_volumes) && return equilibrium_impl!(output, w, x)
    return native_read() do
        try
            equilibrium_impl!(output, w, x)
        catch
            # A failed native solve must not seed the next evaluation.
            w.needs_restore = true
            rethrow()
        end
    end
end

function equilibrium_impl!(output, w::ChemistryWorker, x)
    id = w.id
    n, nb = length(w.reference), size(output, 2)
    np = length(w.molar_volumes)
    if np == 0
        for cell in 1:nb, i in 1:n
            output[i, cell] = x[i, cell]
        end
        output[end, :] .= 0.0
        return output
    end
    if w.needs_restore || !w.reuse_speciation
        checked(RM_StateApply(id, 1), id, "restore baseline")
        w.needs_restore = false
    else
        # Reset all mineral inventories, but retain aqueous activities as
        # nonlinear initial guesses. The complete conserved totals are set
        # below, so these guesses carry no physical or timestep history.
        checked(RM_InitialEquilibriumPhases2Module(id, w.initial_phases), id, "reset minerals")
    end
    for cell in 1:nb
        w.pressure[cell] = x[n + 1, cell]/101325.0
        w.temperature[cell] = x[n + 2, cell] - 273.15
        for i in 1:n
            w.concentrations[cell + (i - 1)*nb] = x[i, cell] + w.reference[i]
        end
    end
    checked(RM_SetPressure(id, w.pressure), id, "pressure")
    checked(RM_SetTemperature(id, w.temperature), id, "temperature")
    checked(RM_SetConcentrations(id, w.concentrations), id, "component totals")
    checked(RM_RunCells(id), id, "equilibrium")
    RM_GetSelectedOutputColumnCount(id) == 2*np || error("Unexpected PHREEQC mineral output layout.")
    checked(RM_GetSelectedOutput(id, w.selected_output), id, "mineral amounts")
    # The native cell is normalized to 1 L at 1 kg/L. Pure-mineral equilibrium
    # is homogeneous in extensive quantities: these amounts are mol/kg carrier,
    # independent of the simulator's carrier density, saturation and volume.
    for cell in 1:nb
        volume = 0.0
        for p in 1:np
            amount = w.selected_output[cell + 2*(p - 1)*nb]
            isfinite(amount) && amount >= 0 || error("Invalid PHREEQC mineral amount.")
            output[n + p, cell] = amount
            w.saturation_indices[p, cell] = w.selected_output[cell + (2*p - 1)*nb]
            volume += amount*w.molar_volumes[p]
        end
        output[end, cell] = volume
        for i in 1:n
            solid = 0.0
            for p in 1:np
                solid += w.stoichiometry[i, p]*output[n + p, cell]
            end
            output[i, cell] = x[i, cell] - solid
        end
    end
    return output
end

function phase_split_jacobian!(w, chemistry, count)
    n, np = length(chemistry.components), length(chemistry.minerals)
    fill!(w.jacobian, 0.0)
    for cell in 1:count, i in 1:n
        w.jacobian[i, i, cell] = 1.0
    end
    # The aqueous map is exactly the identity in the interior of a region
    # without minerals. Near a phase boundary use the one-sided native flash.
    active = [any(p -> w.y[n + p, cell] > 0 || w.saturation_indices[p, cell] >= -1e-7,
        1:np) for cell in 1:count]
    any(active) || return w.jacobian
    nd = n + 1 + Int(chemistry.variable_temperature)
    for j in 1:nd
        copyto!(w.trial, w.x)
        for cell in 1:size(w.x, 2)
            scale = j <= n ? 1e-6 : (j == n + 1 ? 1e5 : 300.0)
            h = chemistry.relative_step*max(abs(w.x[j, cell]), scale)
            w.steps[cell] = h
            w.trial[j, cell] += h
        end
        equilibrium!(w.perturbed, w, w.trial)
        for cell in 1:count, p in 1:np
            if active[cell]
                w.jacobian[n + p, j, cell] = (w.perturbed[n + p, cell] - w.y[n + p, cell])/w.steps[cell]
            end
        end
        if chemistry.difference == :central
            # Central differences where admissible; forward at a zero inventory.
            for cell in 1:size(w.x, 2)
                signed = j <= n && chemistry.signed[j] && chemistry.reference[j] == 0
                inventory = w.x[j, cell] + (j <= n ? chemistry.reference[j] : 0.0)
                central = signed || inventory > w.steps[cell]
                w.trial[j, cell] = w.x[j, cell] + (central ? -w.steps[cell] : w.steps[cell])
            end
            equilibrium!(w.perturbed, w, w.trial)
            for cell in 1:count
                signed = j <= n && chemistry.signed[j] && chemistry.reference[j] == 0
                inventory = w.x[j, cell] + (j <= n ? chemistry.reference[j] : 0.0)
                if active[cell] && (signed || inventory > w.steps[cell])
                    for p in 1:np
                        w.jacobian[n + p, j, cell] = 0.5*(w.jacobian[n + p, j, cell] +
                            (w.y[n + p, cell] - w.perturbed[n + p, cell])/w.steps[cell])
                    end
                end
            end
        end
        for cell in 1:count
            for p in 1:np
                derivative = w.jacobian[n + p, j, cell]
                for i in 1:n
                    w.jacobian[i, j, cell] -= w.stoichiometry[i, p]*derivative
                end
                w.jacobian[end, j, cell] += w.molar_volumes[p]*derivative
            end
        end
    end
    all(isfinite, view(w.jacobian, :, :, 1:count)) || error("Nonfinite PHREEQC phase-split Jacobian.")
    return w.jacobian
end

function close_chemistry!(chemistry)
    lock(chemistry.lifecycle_lock) do
        if !chemistry.closed
            chemistry.closed = true
            while chemistry.active_updates > 0
                wait(chemistry.lifecycle_condition)
            end
            for w in chemistry.workers
                native_write() do
                    checked(RM_Destroy(w.id), w.id, "destroy")
                end
            end
        end
    end
    return nothing
end
