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
    results = Dict{String, Any}()
    # The existing host "fake kernel" benchmark dispatches per cell, whereas
    # production gives each worker a range of cells for all variables.
    if !use_kernel
        results["secondary_serial"] = timer_rows(JB.benchmark_secondary_variables(model, state;
            n = 50, use_kernel = false, fake_kernel = false, verbose = true))
    end
    results["secondary"] = timer_rows(JB.benchmark_secondary_variables(model, state;
        n = 50, use_kernel, fake_kernel = !use_kernel, verbose = true))
    results["assembly"] = timer_rows(JB.benchmark_assembly(sim;
        n = 30, dt = first(case.dt), forces, include_secondary = false))
    results["linear"] = timer_rows(JB.benchmark_linear_operators(sim, config;
        n = 100, n_update = 30, dt = first(case.dt), forces, verbose = true))
    result = Dict("case" => name, "refinement_z" => refinement,
        "cells" => number_of_cells(model.data_domain), "threads" => Threads.nthreads(:default),
        "state" => "initial physical state; existing operation benchmarks; no integration",
        "secondary_uses_kernel" => use_kernel,
        "secondary_note" => "Built-in threaded host loop dispatches per cell; production gives each worker a cell range for all variables",
        "results" => results, "julia" => string(VERSION),
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`))
    open(ARGS[3], "w") do io
        TOML.print(io, result)
    end
end
probe_builtin()
