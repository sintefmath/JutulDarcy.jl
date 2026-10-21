# Paired execution changes with identical Float64 algorithms and fresh simulators.
ENV["PERF_GS_LEVELS"] = "true"
ENV["PERF_GS_RAW"] = "true"
include("benchmark.jl")

function compare_execution()
    name = get(ARGS, 1, "olympus")
    refinement = parse(Int, get(ARGS, 2, "4"))
    @assert get(ENV, "PERF_SINGLE_PRECISION_LINEAR", "false") == "false"
    case = benchmark_case(name, refinement, 3)
    # Typed property tiling now belongs to the library. Earlier four-policy
    # comparisons are recorded at d313001c, against the original library path.
    policies = (("baseline", false), ("gs_raw", true))
    selected = split(get(ENV, "PERF_EXECUTION_POLICIES", "baseline gs_raw"))
    policies = filter(p -> first(p) in selected, policies)
    isempty(policies) && error("Select baseline or gs_raw")
    sources = Dict("jutul" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`))
    function select_policy(levels)
        Jutul.KAPreconditioners.perf_gs_use_levels[] = levels
        Jutul.KAPreconditioners.perf_gs_raw[] = true
    end
    for (label, levels) in policies
        select_policy(levels)
        println("Warming ", label); flush(stdout)
        sim, config = new_simulator(case)
        run_simulation(sim, config, case)
    end
    samples = Any[]
    reference = nothing
    prefix = joinpath(@__DIR__, "results", "execution_$(name)_z$(refinement)_t$(Threads.nthreads(:default))")
    repeats = parse(Int, get(ENV, "PERF_REPEATS", "3"))
    for repeat in 1:repeats
        # Reverse successive rounds to reduce ordering and thermal bias.
        order = isodd(repeat) ? policies : reverse(policies)
        for (label, levels) in order
            select_policy(levels)
            sim, config = new_simulator(case)
            GC.gc()
            measured = @timed run_simulation(sim, config, case)
            stats = report_stats(measured.value.reports)
            state = last(measured.value.states)[:Reservoir]
            current = (pressure = copy(state[:Pressure]), saturations = copy(state[:Saturations]))
            isnothing(reference) && (reference = current)
            @assert all(isfinite, current.pressure) && all(isfinite, current.saturations)
            row = Dict("policy" => label, "repeat" => repeat, "seconds" => measured.time,
                "compile_seconds" => measured.compile_time, "bytes" => measured.bytes,
                "gc_seconds" => measured.gctime, "newtons" => stats.newtons,
                "linear_iterations" => stats.linear_iterations,
                "preconditioner_applications" => stats.linear_precond_iterations,
                "wasted_newtons" => stats.wasted.newtons,
                "pressure_max_absolute" => maximum(abs, current.pressure .- reference.pressure),
                "saturation_max_absolute" => maximum(abs, current.saturations .- reference.saturations),
                "bitwise_equal" => isequal(current, reference),
                "phase_seconds" => Dict(string(k) => v for (k, v) in pairs(stats.time_sum)))
            push!(samples, row)
            println(row); flush(stdout)
            repeat == repeats && serialize(prefix * "_$(label).state", current)
            open(prefix * ".toml", "w") do io
                TOML.print(io, Dict("sources" => sources, "samples" => samples,
                    "float_type" => "Float64", "index_type" => get(ENV, "PERF_INDEX_TYPE", "Int"),
                    "threads" => Threads.nthreads(:default), "case" => name,
                    "refinement_z" => refinement, "report_days" => case.dt ./ si_unit(:day),
                    "pinned_cpus" => get(ENV, "PERF_PIN_CPUS", "")))
            end
            sim = config = measured = nothing
        end
    end
end
compare_execution()
