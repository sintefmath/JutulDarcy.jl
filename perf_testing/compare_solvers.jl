# Compare existing solver configurations in one process with a shared input case.
include("benchmark.jl")

struct ScreeningTimeLimit <: Exception end

function screen_simulator(case; solver_options = NamedTuple())
    sim, config = new_simulator(case; solver_options)
    limit = parse(Float64, get(ENV, "PERF_SCREEN_LIMIT_SECONDS", "Inf"))
    if isfinite(limit)
        started = time_ns()
        config[:post_iteration_hook] = function(converged, report, storage, model, dt, forces, config, iteration)
            (time_ns() - started)/1e9 <= limit || throw(ScreeningTimeLimit())
            return converged
        end
    end
    return sim, config
end

function compare_solvers(; solver_policies = ("default" => NamedTuple(),))
    name = get(ARGS, 1, "olympus")
    refinement = parse(Int, get(ARGS, 2, "2"))
    case = benchmark_case(name, refinement, 3)
    output = Any[]
    solvers = split(get(ENV, "PERF_SOLVERS", "gmres bicgstab"))
    systems = split(get(ENV, "PERF_SYSTEM_SMOOTHERS", "ilu0 ka_ilu0"))
    pressures = split(get(ENV, "PERF_PRESSURE_SMOOTHERS", "spai0 default"))
    repeats = parse(Int, get(ENV, "PERF_REPEATS", "3"))
    label = get(ENV, "PERF_LABEL", "solver_comparison")
    for linear in solvers, system in systems, pressure in pressures, (policy, solver_options) in solver_policies
        ENV["PERF_LINEAR_SOLVER"] = linear
        ENV["PERF_SYSTEM_SMOOTHER"] = system
        ENV["PERF_PRESSURE_SMOOTHER"] = pressure
        println("Warming ", (name, refinement, Threads.nthreads(), linear, system, pressure, policy))
        flush(stdout)
        row = try
            sim, config = screen_simulator(case; solver_options)
            actual_linear_solver = String(config[:linear_solver].solver)
            config[:info_level] = parse(Int, get(ENV, "PERF_SCREEN_INFO", "-1"))
            run_simulation(sim, config, case)
            solver_info = solver_metadata(config)
            samples = Any[]
            for repeat in 1:repeats
                sim, config = screen_simulator(case; solver_options)
                GC.gc()
                measurement = @timed run_simulation(sim, config, case)
                stats = report_stats(measurement.value.reports)
                state = last(measurement.value.states)[:Reservoir]
                push!(samples, Dict("seconds" => measurement.time,
                    "bytes" => measurement.bytes, "gc_seconds" => measurement.gctime,
                    "compile_seconds" => measurement.compile_time,
                    "newtons" => stats.newtons, "linear_iterations" => stats.linear_iterations,
                    "linearizations" => stats.linearizations,
                    "preconditioner_applications" => stats.linear_precond_iterations,
                    "wasted_newtons" => stats.wasted.newtons,
                    "wasted_linear_iterations" => stats.wasted.linear_iterations,
                    "wasted_linearizations" => stats.wasted.linearizations,
                    "wasted_preconditioner_applications" => stats.wasted.linear_precond_iterations,
                    "pressure_sum" => sum(state[:Pressure]),
                    "saturation_norm" => norm(state[:Saturations]),
                    "phase_seconds" => Dict(string(k)=>v for (k,v) in pairs(stats.time_sum))))
                if repeat == repeats && get(ENV, "PERF_SAVE_STATES", "false") == "true"
                    prefix = "$(label)_$(name)_z$(refinement)_t$(Threads.nthreads())_$(policy)_$(system)_$(pressure)_$(linear)"
                    serialize(joinpath(@__DIR__, "results", prefix * ".state"),
                        (pressure = copy(state[:Pressure]), saturations = copy(state[:Saturations])))
                end
                sim = config = measurement = nothing
            end
            Dict("case" => name, "refinement_z" => refinement,
                "threads" => Threads.nthreads(), "linear_solver" => actual_linear_solver,
                "requested_linear_solver" => linear, "solver_options" => repr(solver_options),
                "solver_metadata" => solver_info,
                "system_smoother" => system, "pressure_smoother" => pressure, "policy" => policy,
                "samples" => samples, "median_seconds" => median(s["seconds"] for s in samples))
        catch err
            err isa ScreeningTimeLimit || rethrow()
            Dict("case" => name, "refinement_z" => refinement,
                "threads" => Threads.nthreads(), "linear_solver" => linear,
                "system_smoother" => system, "pressure_smoother" => pressure, "policy" => policy,
                "status" => "screening_time_limit", "limit_seconds" =>
                    parse(Float64, ENV["PERF_SCREEN_LIMIT_SECONDS"]))
        end
        push!(output, row)
        println(row)
        flush(stdout)
        open(joinpath(@__DIR__, "results", "$(label)_$(name)_z$(refinement)_t$(Threads.nthreads()).toml"), "w") do io
            TOML.print(io, Dict("results" => output,
                "report_days" => case.dt ./ si_unit(:day),
                "julia" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
                "screen_limit_seconds" => get(ENV, "PERF_SCREEN_LIMIT_SECONDS", "Inf"),
                "single_precision_linear_solver" => get(ENV, "PERF_SINGLE_PRECISION_LINEAR", "false"),
                "index_type" => get(ENV, "PERF_INDEX_TYPE", "Int"),
                "relaxation" => "false",
                "pinned_cpus" => get(ENV, "PERF_PIN_CPUS", ""),
                "gs_levels" => get(ENV, "PERF_GS_LEVELS", "false"),
                "gs_global_chunks" => get(ENV, "PERF_GS_GLOBAL_CHUNKS", "false"),
                "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
                "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`)))
        end
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    compare_solvers()
end
