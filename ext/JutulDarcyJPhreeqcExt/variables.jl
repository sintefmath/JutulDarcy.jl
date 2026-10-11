function Jutul.get_dependencies(c::PhreeqcEquilibrium, model)
    deps = (:TracerConcentrations, :Pressure, :Saturations)
    return c.variable_temperature ? (deps..., :Temperature) : deps
end

Jutul.get_dependencies(::PhreeqcSplitConcentrations, model) = (:PhreeqcPhaseSplit, )
Jutul.get_dependencies(::PhreeqcMineralVolume, model) = (:PhreeqcPhaseSplit, )

Jutul.@jutul_secondary function update_split_concentrations!(out, v::PhreeqcSplitConcentrations,
        model, PhreeqcPhaseSplit, ix)
    for cell in ix, (i, row) in enumerate(v.rows)
        out[i, cell] = PhreeqcPhaseSplit[row, cell]
    end
    return out
end

Jutul.@jutul_secondary function update_mineral_volume!(out, v::PhreeqcMineralVolume,
        model, PhreeqcPhaseSplit, ix)
    for cell in ix
        out[cell] = PhreeqcPhaseSplit[v.row, cell]
    end
    return out
end

function Jutul.get_dependencies(v::PhreeqcPoreVolume, model)
    return (v.base_volume, :PhreeqcReferenceMineralVolume, :PhreeqcMineralVolume,
        :PhaseMassDensities, :Saturations)
end

function Jutul.update_secondary_variable!(out, v::PhreeqcPoreVolume, model, state, ix)
    base = state[v.base_volume]
    reference = state.PhreeqcReferenceMineralVolume
    for cell in ix
        density = state.PhaseMassDensities[v.phase_index, cell]
        saturation = state.Saturations[v.phase_index, cell]
        # V = V_base + V_s,ref - V_s, with V_s = rho*S*V*nu_s.
        # Solve this local relation explicitly; no circular dependencies,
        # lagged porosity, timestep history or numerical volume derivative.
        out[cell] = (base[cell] + reference[cell])/
            (1 + density*saturation*state.PhreeqcMineralVolume[cell])
        ForwardDiff.value(out[cell]) > 0 || error("Nonpositive reactive pore volume in cell $cell.")
    end
    return out
end

Jutul.get_dependencies(::PhreeqcImmobileConcentrations, model) =
    (:TracerConcentrations, :PhreeqcAqueousConcentrations)

Jutul.@jutul_secondary function update_immobile!(out, c::PhreeqcImmobileConcentrations,
        model, TracerConcentrations, PhreeqcAqueousConcentrations, ix)
    for cell in ix, i in 1:c.ncomponents
        out[i, cell] = TracerConcentrations[i, cell] - PhreeqcAqueousConcentrations[i, cell]
    end
    return out
end

@inline function chemistry_input(state, chemistry, j, cell)
    n = length(chemistry.components)
    if j <= n
        return state.TracerConcentrations[j, cell]
    elseif j == n + 1
        return state.Pressure[cell]
    else
        return chemistry.variable_temperature ? state.Temperature[cell] : chemistry.temperature
    end
end

function validate_input(x, chemistry, cell)
    n = length(chemistry.components)
    all(isfinite, x) || error("Nonfinite PHREEQC input in cell $cell.")
    for i in 1:n
        if !chemistry.signed[i] && x[i] < 0.0
            error("Negative total inventory of $(chemistry.components[i]) in cell $cell.")
        end
        if chemistry.reference[i] > 0 && x[i] + chemistry.reference[i] <= 0
            error("PHREEQC requires positive total $(chemistry.components[i]) inventory in cell $cell.")
        end
    end
    x[n + 1] > 0 || error("PHREEQC requires positive absolute pressure in cell $cell.")
    x[n + 2] > 0 || error("PHREEQC requires positive absolute temperature in cell $cell.")
end

function write_split!(out, chemistry, state, cell, derivatives)
    nd = size(chemistry.cache_input, 1)
    for i in axes(out, 1)
        v = chemistry.cache_output[i, cell]
        if derivatives
            # Lift the native Float64 result through the full local Jacobian.
            # This propagates cross-component, pressure and temperature
            # partials. Density enters the volume closure directly, rather
            # than through the normalized pure-mineral equilibrium map.
            dv = zero(out[i, cell])
            for j in 1:nd
                x = chemistry_input(state, chemistry, j, cell)
                dv += chemistry.cache_jacobian[i, j, cell]*(x - ForwardDiff.value(x))
            end
            out[i, cell] = v + dv
        else
            out[i, cell] = v
        end
    end
end

function update_batch!(out, chemistry, state, cells, derivatives)
    worker = chemistry_worker(chemistry)
    no, nd = size(chemistry.cache_output, 1), size(chemistry.cache_input, 1)
    dirty = Int[]
    for cell in cells
        changed = !chemistry.cache_valid[cell]
        for j in 1:nd
            x = Float64(ForwardDiff.value(chemistry_input(state, chemistry, j, cell)))
            changed |= x != chemistry.cache_input[j, cell]
        end
        if changed || (derivatives && !chemistry.cache_derivatives[cell])
            push!(dirty, cell)
        else
            write_split!(out, chemistry, state, cell, derivatives)
        end
    end
    count = length(dirty)
    count == 0 && return
    for (local_cell, cell) in enumerate(dirty)
        for j in 1:nd
            worker.x[j, local_cell] = Float64(ForwardDiff.value(chemistry_input(state, chemistry, j, cell)))
        end
        validate_input(view(worker.x, :, local_cell), chemistry, cell)
    end
    # Pad the native batch with the last valid state. The padding is never
    # cached or written to Jutul; it avoids resizing any native-side arrays.
    for local_cell in (count + 1):chemistry.batch_size, j in 1:nd
        worker.x[j, local_cell] = worker.x[j, count]
    end
    try
        equilibrium!(worker.y, worker, worker.x)
        derivatives && phase_split_jacobian!(worker, chemistry, count)
    catch err
        error("PHREEQC update for cells $(first(dirty)):$(last(dirty)) failed: $(sprint(showerror, err))")
    end
    for (local_cell, cell) in enumerate(dirty)
        for j in 1:nd
            chemistry.cache_input[j, cell] = worker.x[j, local_cell]
            if derivatives
                for i in 1:no
                    chemistry.cache_jacobian[i, j, cell] = worker.jacobian[i, j, local_cell]
                end
            end
        end
        for i in 1:no
            chemistry.cache_output[i, cell] = worker.y[i, local_cell]
        end
        chemistry.cache_derivatives[cell] = derivatives
        chemistry.cache_valid[cell] = true
        write_split!(out, chemistry, state, cell, derivatives)
    end
end

function Jutul.update_secondary_variable!(out, chemistry::PhreeqcEquilibrium, model, state, ix)
    lock(chemistry.lifecycle_lock) do
        chemistry.closed && error("The PHREEQC runtime has been closed. Set up a new model before using chemistry.")
        chemistry.active_updates += 1
    end
    try
        size(out) == size(chemistry.cache_output) || error("PHREEQC phase-split output dimensions do not match the model.")
        cells = collect(ix)
        isempty(cells) && return out
        for cell in cells
            ForwardDiff.value(state.Saturations[chemistry.phase_index, cell]) > 1e-12 ||
                error("PHREEQC requires a present aqueous carrier in cell $cell; dry cells require a different inventory formulation.")
        end
        derivatives = Jutul.number_of_partials(eltype(out)) > 0
        # Global batch locks synchronize overlapping model-copy evaluations;
        # Jutul's disjoint cell partitions can evaluate concurrently.
        groups = Dict{Int, Vector{Int}}()
        for cell in cells
            bucket = cld(cell, chemistry.batch_size)
            push!(get!(() -> Int[], groups, bucket), cell)
        end
        buckets = sort!(collect(keys(groups)))
        nbatches = length(buckets)
        ntasks = min(length(chemistry.workers), nbatches)
        function update_part(task_index)
            with_chemistry_worker(chemistry) do worker
                for batch in task_index:ntasks:nbatches
                    bucket = buckets[batch]
                    lock(chemistry.cache_locks[bucket]) do
                        update_batch!(out, chemistry, state, groups[bucket], derivatives)
                    end
                end
            end
        end
        # Jutul already threads partial cell ranges, potentially with Polyester.
        # Avoid nesting Julia tasks inside that loop; full-range calls manage
        # their own parallelism so even small flow models thread chemistry.
        if ntasks == 1 || length(cells) != size(out, 2)
            ntasks = 1
            update_part(1)
        else
            # Ownership follows tasks, not threadid(): safe under task migration
            # and when Jutul also partitions secondary updates across threads.
            @sync for task_index in 1:ntasks
                Threads.@spawn update_part(task_index)
            end
        end
    finally
        lock(chemistry.lifecycle_lock) do
            chemistry.active_updates -= 1
            notify(chemistry.lifecycle_condition; all = true)
        end
    end
    return out
end
