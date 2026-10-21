# Usage: julia --project=perf_testing --threads=10,0 perf_testing/benchmark.jl olympus 1 3 3 baseline
using Jutul, JutulDarcy, LinearAlgebra, Statistics, TOML, Serialization, Profile

include("refine_layers.jl")

BLAS.set_num_threads(1)

function benchmark_case(name, refinement, steps)
    dataset = Dict("norne" => "NORNE_NOHYST", "olympus" => "OLYMPUS_1")[name]
    path = JutulDarcy.GeoEnergyIO.test_input_file_path(dataset, "$dataset.DATA")
    data = JutulDarcy.GeoEnergyIO.parse_data_file(path)
    if refinement > 1
        data = refine_layers(data, refinement)
    end
    case = setup_case_from_data_file(data; backend = :csr, block_backend = true)
    short = case[1:min(steps, length(case.dt))]
    # Olympus uses annual report steps: shorten the physical interval as well.
    max_days = parse(Float64, get(ENV, "PERF_MAX_DAYS", "1"))
    short.dt .= min.(short.dt, max_days * si_unit(:day))
    return short
end

function new_simulator(case)
    smoother = Symbol(get(ENV, "PERF_PRESSURE_SMOOTHER", "spai0"))
    linear_solver_arg = smoother == :default ? NamedTuple() : (amg_arg = (smoother_type = smoother,),)
    system_smoother = Symbol(get(ENV, "PERF_SYSTEM_SMOOTHER", "default"))
    if system_smoother != :default
        linear_solver_arg = merge(linear_solver_arg, (smoother_type = system_smoother,))
    end
    return setup_reservoir_simulator(case;
        mode = :ka, linear_solver = :gmres, info_level = -1,
        linear_solver_arg, output_substates = false)
end

function run_simulation(sim, config, case)
    result = simulate!(sim, case.dt; forces = case.forces, state0 = case.state0, config)
    @assert length(result.states) == length(case.dt) "Simulation did not complete"
    for report in result.reports
        @assert last(report[:ministeps])[:success] "Report interval did not complete"
    end
    return result
end

function main()
    name = get(ARGS, 1, "olympus")
    refinement = parse(Int, get(ARGS, 2, "1"))
    steps = parse(Int, get(ARGS, 3, "3"))
    repeats = parse(Int, get(ARGS, 4, "3"))
    label = get(ARGS, 5, "baseline")
    outdir = get(ENV, "PERF_OUTPUT", joinpath(@__DIR__, "results"))
    mkpath(outdir)
    println("Preparing ", (name, refinement, steps, Threads.nthreads(:default)))
    flush(stdout)
    case = benchmark_case(name, refinement, steps)
    nc = number_of_cells(reservoir_model(case.model).data_domain)
    println("Warming up ", nc, " cells; report days = ", case.dt ./ si_unit(:day))
    flush(stdout)
    sim, config = new_simulator(case)
    run_simulation(sim, config, case)
    sim = config = nothing
    samples = Any[]
    prefix = joinpath(outdir, "$(label)_$(name)_z$(refinement)_t$(Threads.nthreads(:default))")
    for repeat in 1:repeats
        # Fresh hierarchy for every sample: warm code, same initial physical state.
        sim, config = new_simulator(case)
        GC.gc()
        measurement = @timed run_simulation(sim, config, case)
        stats = report_stats(measurement.value.reports)
        state = last(measurement.value.states)[:Reservoir]
        pressure = state[:Pressure]
        saturations = state[:Saturations]
        sample = Dict{String, Any}(
            "repeat" => repeat, "seconds" => measurement.time,
            "gc_seconds" => measurement.gctime, "bytes" => measurement.bytes,
            "newtons" => stats.newtons, "linear_iterations" => stats.linear_iterations,
            "ministeps" => stats.ministeps,
            "wasted_newtons" => stats.wasted.newtons,
            "wasted_linear_iterations" => stats.wasted.linear_iterations,
            "pressure_sum" => sum(pressure), "pressure_norm" => norm(pressure),
            "saturation_norm" => norm(saturations),
            "phase_seconds" => Dict(string(k) => v for (k, v) in pairs(stats.time_sum)))
        push!(samples, sample)
        println(sample)
        flush(stdout)
        if repeat == repeats
            serialize(prefix * ".state", (pressure = copy(pressure), saturations = copy(saturations)))
        end
        sim = config = measurement = nothing
    end
    result = Dict(
        "report_days" => case.dt ./ si_unit(:day),
        "case" => name, "refinement_z" => refinement, "report_steps" => steps,
        "cells" => nc, "threads" => Threads.nthreads(:default),
        "pressure_smoother" => get(ENV, "PERF_PRESSURE_SMOOTHER", "spai0"),
        "system_smoother" => get(ENV, "PERF_SYSTEM_SMOOTHER", "default"),
        "julia" => string(VERSION), "label" => label,
        "package_versions" => Dict("Jutul" => string(pkgversion(Jutul)),
            "JutulDarcy" => string(pkgversion(JutulDarcy)),
            "Krylov" => string(pkgversion(Jutul.Krylov)),
            "KernelAbstractions" => string(pkgversion(Jutul.KernelExecution.KernelAbstractions)),
            "Polyester" => string(pkgversion(Jutul.Polyester))),
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`),
        "samples" => samples,
        "median_seconds" => median(s["seconds"] for s in samples))
    open(prefix * ".toml", "w") do io
        TOML.print(io, result)
    end
    if get(ENV, "PERF_PROFILE", "false") == "true"
        sim, config = new_simulator(case)
        GC.gc()
        Profile.clear()
        @profile run_simulation(sim, config, case)
        open(prefix * ".profile.txt", "w") do io
            Profile.print(io; format = :flat, sortedby = :count, mincount = 20, C = true)
        end
        serialize(prefix * ".profile", Profile.retrieve())
    end
    println("Saved ", prefix, ".toml; median = ", result["median_seconds"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
