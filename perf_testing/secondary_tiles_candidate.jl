# Process-local host experiment, included into Jutul. Keep worker limits and
# variable order; use the existing minbatch as the cache-tile size by default.
function update_secondary_variables_state!(state, model::SimulationModel, vars::NamedTuple)
    ctx = model.context
    eligible = ctx isa KernelAbstractionsContext && KernelExecution.is_cpu_backend(ctx) &&
        !ctx.use_kernels_for_secondary && Threads.nthreads(:default) > 1 &&
        nthreads(ctx) > 1 && thread_type(ctx) == :batch && !isempty(vars)
    cells = eligible ? number_of_entities(model, first(values(vars))) : 0
    if !eligible || cells < 16_384 ||
            !all(number_of_entities(model, var) == cells for var in values(vars))
        return invoke(update_secondary_variables_state!, Tuple{Any, Any, Any}, state, model, vars)
    end
    workers = clamp(cells ÷ minbatch(ctx), 1, min(nthreads(ctx), Threads.nthreads(:default)))
    tile_size = parse(Int, get(ENV, "PERF_SECONDARY_TILE", string(minbatch(ctx))))
    tile_size > 0 || throw(ArgumentError("PERF_SECONDARY_TILE must be positive"))
    tiles = max(workers, cld(cells, tile_size))
    definitions = map((symbol, variable) -> (Val(symbol), variable), keys(vars), values(vars))
    perf_update_secondary_tiles!(state, model, definitions, workers, tiles, thread_type(ctx))
    return state
end

@inline function perf_update_secondary_tile!(state, model, ::Val{K}, variable, tile, tiles) where K
    target = state[K]
    indices = entity_eachindex(target, tile, tiles)
    update_secondary_variable!(target, variable, model, state, indices)
    return nothing
end

function perf_update_secondary_tiles!(state, model, definitions::Tuple, workers, tiles, threading)
    function update_worker(worker)
        for tile in load_balanced_interval(worker, tiles, workers)
            foreach(definitions) do (symbol, variable)
                perf_update_secondary_tile!(state, model, symbol, variable, tile, tiles)
            end
        end
    end
    threaded_loop(update_worker, workers, threading)
    return nothing
end
