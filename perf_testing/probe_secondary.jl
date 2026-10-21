# Property-update scaling on the same initial physical state; no timestepping.
include("benchmark.jl")

function tiled_secondary!(state, model, variables, tile_size)
    context = model.context
    first_variable = last(first(pairs(variables)))
    cells = Jutul.number_of_entities(model, first_variable)
    workers = clamp(cells ÷ Jutul.minbatch(context), 1,
        min(Jutul.nthreads(context), Threads.nthreads(:default)))
    tiles = max(workers, cld(cells, tile_size))
    function update_worker(worker)
        for tile in Jutul.load_balanced_interval(worker, tiles, workers)
            for (symbol, variable) in pairs(variables)
                values = state[symbol]
                indices = Jutul.entity_eachindex(values, tile, tiles)
                Jutul.update_secondary_variable!(values, variable, model, state, indices)
            end
        end
    end
    Jutul.threaded_loop(update_worker, workers, Jutul.thread_type(context))
    return nothing
end

function secondary_measure(f::F) where F
    f()
    samples = [@timed(for _ in 1:5; f(); end) for _ in 1:9]
    return Dict("seconds" => median(s.time for s in samples)/5,
        "bytes" => minimum(s.bytes for s in samples)/5,
        "compile_seconds" => sum(s.compile_time for s in samples))
end

function probe_secondary()
    name, refinement = ARGS[1], parse(Int, ARGS[2])
    case = benchmark_case(name, refinement, 3)
    sim, _ = new_simulator(case)
    storage, model = Jutul.submodel_evaluation_pair(sim.storage, sim.model, :Reservoir)
    variables = storage.variable_definitions.secondary_variables
    original() = Jutul.update_secondary_variables!(storage, model)
    original()
    state = Jutul.evaluation_state(storage)
    first_variable = last(first(pairs(variables)))
    cells = Jutul.number_of_entities(model, first_variable)
    @assert all(Jutul.number_of_entities(model, variable) == cells for (_, variable) in pairs(variables))
    reference = Dict(symbol => copy(state[symbol]) for (symbol, _) in pairs(variables))
    @assert all(p -> isfinite(p) && Jutul.value(p) > 0, state[:Pressure])
    rows = Any[merge(Dict("tile" => "original"), secondary_measure(original))]
    for size in (256, 512, 1024, 2048, 4096, 8192)
        candidate = let tile_size = size
            () -> tiled_secondary!(state, model, variables, tile_size)
        end
        candidate()
        @assert all(isequal(state[symbol], expected) for (symbol, expected) in reference)
        row = merge(Dict("tile" => string(size), "bitwise_equal" => true), secondary_measure(candidate))
        push!(rows, row)
        println(row)
        flush(stdout)
    end
    result = Dict("case" => name, "refinement_z" => refinement, "cells" => cells,
        "threads" => Threads.nthreads(:default), "context_workers" => Jutul.nthreads(model.context),
        "state" => "initial condition, no timestep integration", "variables" => string.(keys(reference)),
        "rows" => rows, "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`))
    open(ARGS[3], "w") do io
        TOML.print(io, result)
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    probe_secondary()
end
