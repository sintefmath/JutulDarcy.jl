# Compare existing solver configurations in one process with a shared input case.
include("benchmark.jl")

function compare_solvers()
    name = get(ARGS, 1, "olympus")
    refinement = parse(Int, get(ARGS, 2, "2"))
    case = benchmark_case(name, refinement, 3)
    output = Any[]
    solvers = split(get(ENV, "PERF_SOLVERS", "gmres bicgstab"))
    systems = split(get(ENV, "PERF_SYSTEM_SMOOTHERS", "ilu0 ka_ilu0"))
    pressures = split(get(ENV, "PERF_PRESSURE_SMOOTHERS", "spai0 default"))
    repeats = parse(Int, get(ENV, "PERF_REPEATS", "3"))
    label = get(ENV, "PERF_LABEL", "solver_comparison")
    for linear in solvers, system in systems, pressure in pressures
        ENV["PERF_LINEAR_SOLVER"] = linear
        ENV["PERF_SYSTEM_SMOOTHER"] = system
        ENV["PERF_PRESSURE_SMOOTHER"] = pressure
        println("Warming ", (name, refinement, Threads.nthreads(), linear, system, pressure))
        flush(stdout)
        sim, config = new_simulator(case)
        run_simulation(sim, config, case)
        samples = Any[]
        for _ in 1:repeats
            sim, config = new_simulator(case)
            GC.gc()
            measurement = @timed run_simulation(sim, config, case)
            stats = report_stats(measurement.value.reports)
            state = last(measurement.value.states)[:Reservoir]
            push!(samples, Dict("seconds" => measurement.time,
                "bytes" => measurement.bytes, "gc_seconds" => measurement.gctime,
                "compile_seconds" => measurement.compile_time,
                "newtons" => stats.newtons, "linear_iterations" => stats.linear_iterations,
                "wasted_newtons" => stats.wasted.newtons,
                "wasted_linear_iterations" => stats.wasted.linear_iterations,
                "pressure_sum" => sum(state[:Pressure]),
                "saturation_norm" => norm(state[:Saturations]),
                "phase_seconds" => Dict(string(k)=>v for (k,v) in pairs(stats.time_sum))))
            sim = config = measurement = nothing
        end
        row = Dict("case" => name, "refinement_z" => refinement,
            "threads" => Threads.nthreads(), "linear_solver" => linear,
            "system_smoother" => system, "pressure_smoother" => pressure,
            "samples" => samples, "median_seconds" => median(s["seconds"] for s in samples))
        push!(output, row)
        println(row)
        flush(stdout)
        open(joinpath(@__DIR__, "results", "$(label)_$(name)_z$(refinement)_t$(Threads.nthreads()).toml"), "w") do io
            TOML.print(io, Dict("results" => output,
                "report_days" => case.dt ./ si_unit(:day),
                "julia" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
                "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
                "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`)))
        end
    end
end
compare_solvers()
