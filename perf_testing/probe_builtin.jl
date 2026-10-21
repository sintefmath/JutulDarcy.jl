# Existing Jutul operation benchmarks on a fixed initial physical state.
include("benchmark.jl")
const JB = Jutul.Benchmarking
const TO = JB.TimerOutputs

function timer_rows(timer, prefix = "")
    rows = Any[]
    for (name, child) in sort!(collect(timer.inner_timers); by = first)
        path = isempty(prefix) ? name : "$prefix/$name"
        calls = TO.ncalls(child)
        push!(rows, Dict("term" => path, "calls" => calls,
            "seconds" => TO.time(child)/1e9, "bytes" => TO.allocated(child),
            "seconds_per_call" => TO.time(child)/1e9/max(calls, 1)))
        append!(rows, timer_rows(child, path))
    end
    return rows
end

function probe_builtin()
    name, refinement = ARGS[1], parse(Int, ARGS[2])
    case = benchmark_case(name, refinement, 3)
    cells = number_of_cells(reservoir_model(case.model).data_domain)
    sim, config = new_simulator(case)
    forces = case.forces isa AbstractVector ? first(case.forces) : case.forces
    Jutul.update_secondary_variables!(sim.storage, sim.model, true)
    Jutul.initialize_before_first_timestep!(sim, first(case.dt); forces, config)
    Jutul.update_before_step!(sim, first(case.dt), forces; time = 0.0,
        recorder = Jutul.progress_recorder(sim))
    storage, model = Jutul.submodel_evaluation_pair(sim.storage, sim.model, :Reservoir)
    state = Jutul.evaluation_state(storage)
    ctx = model.context
    use_kernel = ctx isa Jutul.KernelAbstractionsContext && ctx.use_kernels_for_secondary
    sections = split(get(ENV, "PERF_BUILTIN_SECTIONS", "secondary assembly linear"))
    results = Dict{String, Any}()
    result = Dict("case" => name, "refinement_z" => refinement,
        "cells" => cells, "threads" => Threads.nthreads(:default),
        "state" => "initial physical state; existing operation benchmarks; no integration",
        "secondary_uses_kernel" => use_kernel,
        "sections" => sections, "reduce_memory" => ctx.reduce_memory,
        "secondary_note" => "Built-in threaded host loop dispatches per cell; production gives each worker a cell range for all variables",
        "results" => results, "julia" => string(VERSION), "complete" => false,
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`))
    function save_progress()
        open(ARGS[3], "w") do io
            TOML.print(io, result)
        end
    end
    # The existing host "fake kernel" benchmark dispatches per cell, whereas
    # production gives each worker a range of cells for all variables.
    if "secondary" in sections && !use_kernel
        results["secondary_serial"] = timer_rows(JB.benchmark_secondary_variables(model, state;
            n = 50, use_kernel = false, fake_kernel = false, verbose = true))
        save_progress()
    end
    if "secondary" in sections
        results["secondary"] = timer_rows(JB.benchmark_secondary_variables(model, state;
            n = 50, use_kernel, fake_kernel = !use_kernel, verbose = true))
        save_progress()
    end
    if "assembly" in sections
        results["assembly"] = timer_rows(JB.benchmark_assembly(sim;
            n = 30, dt = first(case.dt), forces, include_secondary = false))
        save_progress()
    end
    if "linear" in sections
        results["linear"] = timer_rows(JB.benchmark_linear_operators(sim, config;
            n = 100, n_update = 30, dt = first(case.dt), forces, verbose = true))
    end
    result["complete"] = true
    save_progress()
end
probe_builtin()
