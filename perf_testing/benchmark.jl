# Usage: julia --project=perf_testing --threads=10,0 perf_testing/benchmark.jl olympus 1 3 3 baseline
using Jutul, JutulDarcy, LinearAlgebra, Statistics, TOML, Serialization, Profile, SHA

include("refine_layers.jl")
include("pin_threads.jl")
if parse(Bool, get(ENV, "PERF_GS_LEVELS", "false"))
    Base.include(Jutul.KAPreconditioners, joinpath(@__DIR__, "gs_levels_candidate.jl"))
    Jutul.KAPreconditioners.perf_gs_use_levels[] = true
    Jutul.KAPreconditioners.perf_gs_global_chunks[] = parse(Bool, get(ENV, "PERF_GS_GLOBAL_CHUNKS", "false"))
end

BLAS.set_num_threads(parse(Int, get(ENV, "PERF_BLAS_THREADS", "1")))

function benchmark_case(name, refinement, steps)
    dataset = Dict("norne" => "NORNE_NOHYST", "olympus" => "OLYMPUS_1")[name]
    path = JutulDarcy.GeoEnergyIO.test_input_file_path(dataset, "$dataset.DATA")
    data = JutulDarcy.GeoEnergyIO.parse_data_file(path)
    if refinement > 1
        data = refine_layers(data, refinement)
    end
    case = setup_case_from_data_file(data; backend = :csr, block_backend = true)
    fraction = parse(Float64, get(ENV, "PERF_FRACTION", "0"))
    0 <= fraction <= 1 || throw(ArgumentError("PERF_FRACTION must be in [0, 1]"))
    if fraction > 0
        steps = ceil(Int, fraction * length(case.dt))
    end
    short = case[1:min(steps, length(case.dt))]
    # Olympus uses annual report steps: shorten the physical interval as well.
    max_days = parse(Float64, get(ENV, "PERF_MAX_DAYS", fraction > 0 ? "Inf" : "1"))
    short.dt .= min.(short.dt, max_days * si_unit(:day))
    return short
end

function new_simulator(case; solver_options = NamedTuple())
    smoother = Symbol(get(ENV, "PERF_PRESSURE_SMOOTHER", "spai0"))
    linear_solver_arg = smoother == :default ? NamedTuple() : (amg_arg = (smoother_type = smoother,),)
    system_smoother = Symbol(get(ENV, "PERF_SYSTEM_SMOOTHER", "default"))
    if system_smoother != :default
        linear_solver_arg = merge(linear_solver_arg, (smoother_type = system_smoother,))
    end
    linear_solver_arg = merge(linear_solver_arg, solver_options)
    linear_solver = Symbol(get(ENV, "PERF_LINEAR_SOLVER", "gmres"))
    mode = Symbol(get(ENV, "PERF_MODE", "ka"))
    precond = Symbol(get(ENV, "PERF_PRECONDITIONER", "cpr"))
    precond in (:cpr, :cprw) || throw(ArgumentError("Select cpr or cprw"))
    single_precision_linear_solver = parse(Bool, get(ENV, "PERF_SINGLE_PRECISION_LINEAR", "false"))
    index_type = Dict("Int" => Int, "Int32" => Int32)[get(ENV, "PERF_INDEX_TYPE", "Int")]
    reduce_memory = parse(Bool, get(ENV, "PERF_REDUCE_MEMORY", "true"))
    return setup_reservoir_simulator(case;
        mode, linear_solver, precond, info_level = -1,
        linear_solver_arg, single_precision_linear_solver, index_type, reduce_memory,
        relaxation = false, output_substates = false)
end

function run_simulation(sim, config, case)
    result = simulate!(sim, case.dt; forces = case.forces, state0 = case.state0, config)
    @assert length(result.states) == length(case.dt) "Simulation did not complete"
    for report in result.reports
        @assert last(report[:ministeps])[:success] "Report interval did not complete"
    end
    return result
end

function solver_metadata(config)
    preconditioner = config[:linear_solver].preconditioner
    pressure = preconditioner.pressure_precond
    system = preconditioner.system_precond
    metadata = Dict{String,Any}(
        "preconditioner_variant" => String(preconditioner.variant),
        "linear_solver" => String(config[:linear_solver].solver),
        "pressure_preconditioner" => string(typeof(pressure)),
        "system_preconditioner" => string(typeof(system)))
    if hasproperty(system, :factor)
        metadata["system_factor_type"] = string(typeof(system.factor))
    end
    if pressure isa Jutul.AMGPreconditioner
        metadata["pressure_smoother"] = repr(pressure.options.smoother)
        metadata["amg_level_rows"] = [size(level.A, 1) for level in pressure.factor.levels]
        metadata["amg_level_nonzeros"] = [length(level.A.nzval) for level in pressure.factor.levels]
    end
    return metadata
end

function main(; solver_options = NamedTuple())
    name = get(ARGS, 1, "olympus")
    refinement = parse(Int, get(ARGS, 2, "1"))
    steps = parse(Int, get(ARGS, 3, "3"))
    repeats = parse(Int, get(ARGS, 4, "3"))
    label = get(ARGS, 5, "baseline")
    source_state = Dict(string(nameof(package)) => Dict(
        "commit" => readchomp(`git -C $(pkgdir(package)) rev-parse HEAD`),
        "source_diff_sha256" => bytes2hex(sha256(read(`git -C $(pkgdir(package)) diff HEAD -- src`))))
        for package in (Jutul, JutulDarcy))
    outdir = get(ENV, "PERF_OUTPUT", joinpath(@__DIR__, "results"))
    mkpath(outdir)
    println("Preparing ", (name, refinement, steps, Threads.nthreads(:default)))
    flush(stdout)
    case = benchmark_case(name, refinement, steps)
    nc = number_of_cells(reservoir_model(case.model).data_domain)
    println("Warming up ", nc, " cells; report days = ", case.dt ./ si_unit(:day))
    flush(stdout)
    sim, config = new_simulator(case; solver_options)
    warm_case = case
    if haskey(ENV, "PERF_WARMUP_MAX_DAYS")
        warm_steps = parse(Int, get(ENV, "PERF_WARMUP_STEPS", "3"))
        warm_case = case[1:min(warm_steps, length(case.dt))]
        warm_case.dt .= min.(warm_case.dt, parse(Float64, ENV["PERF_WARMUP_MAX_DAYS"])*si_unit(:day))
    end
    # Show physical progress for long warm-ups; timed configurations stay quiet.
    if parse(Float64, get(ENV, "PERF_FRACTION", "0")) > 0
        config[:info_level] = parse(Int, get(ENV, "PERF_WARMUP_INFO", "1"))
    end
    run_simulation(sim, config, warm_case)
    println("Warm-up complete; starting measured runs")
    flush(stdout)
    warmup_solver_info = solver_metadata(config)
    solver_info = warmup_solver_info
    sim = config = nothing
    samples = Any[]
    prefix = joinpath(outdir, "$(label)_$(name)_z$(refinement)_t$(Threads.nthreads(:default))")
    for repeat in 1:repeats
        # Fresh hierarchy for every sample: warm code, same initial physical state.
        sim, config = new_simulator(case; solver_options)
        config[:info_level] = parse(Int, get(ENV, "PERF_RUN_INFO", "-1"))
        GC.gc()
        measurement = @timed run_simulation(sim, config, case)
        stats = report_stats(measurement.value.reports)
        state = last(measurement.value.states)[:Reservoir]
        pressure = state[:Pressure]
        saturations = state[:Saturations]
        sample = Dict{String, Any}(
            "repeat" => repeat, "seconds" => measurement.time,
            "gc_seconds" => measurement.gctime, "bytes" => measurement.bytes,
            "compile_seconds" => measurement.compile_time,
            "recompile_seconds" => measurement.recompile_time,
            "newtons" => stats.newtons, "linear_iterations" => stats.linear_iterations,
            "linearizations" => stats.linearizations,
            "preconditioner_applications" => stats.linear_precond_iterations,
            "ministeps" => stats.ministeps,
            "wasted_newtons" => stats.wasted.newtons,
            "wasted_linear_iterations" => stats.wasted.linear_iterations,
            "wasted_linearizations" => stats.wasted.linearizations,
            "wasted_preconditioner_applications" => stats.wasted.linear_precond_iterations,
            "pressure_sum" => sum(pressure), "pressure_norm" => norm(pressure),
            "saturation_norm" => norm(saturations),
            "phase_seconds" => Dict(string(k) => v for (k, v) in pairs(stats.time_sum)))
        push!(samples, sample)
        println(sample)
        flush(stdout)
        if repeat == repeats
            solver_info = solver_metadata(config)
            serialize(prefix * ".state", (pressure = copy(pressure), saturations = copy(saturations)))
            mesh = Jutul.physical_representation(reservoir_model(case.model).data_domain)
            if hasproperty(mesh, :cell_map)
                serialize(prefix * ".cells", copy(mesh.cell_map))
            end
            if get(ENV, "PERF_SAVE_FIXTURE", "false") == "true"
                serialize(prefix * ".fixture", (
                    preconditioner = config[:linear_solver].preconditioner,
                    linear_system = sim.storage.LinearizedSystem))
            end
        end
        sim = config = measurement = nothing
    end
    result = Dict(
        "report_days" => case.dt ./ si_unit(:day),
        "warmup_report_days" => warm_case.dt ./ si_unit(:day),
        "run_info_level" => get(ENV, "PERF_RUN_INFO", "-1"),
        "case" => name, "refinement_z" => refinement, "report_steps" => length(case.dt),
        "schedule_fraction" => parse(Float64, get(ENV, "PERF_FRACTION", "0")),
        "cells" => nc, "threads" => Threads.nthreads(:default),
        "phase_names" => collect(string.(JutulDarcy.phase_names(reservoir_model(case.model).system))),
        "blas_threads" => BLAS.get_num_threads(),
        "cpu_set" => get(ENV, "PERF_CPUSET", "unrestricted"),
        "omp_threads" => get(ENV, "OMP_NUM_THREADS", "runtime default"),
        "omp_wait_policy" => get(ENV, "OMP_WAIT_POLICY", "runtime default"),
        "solver_metadata" => solver_info,
        "warmup_solver_metadata" => warmup_solver_info,
        "solver_options" => repr(solver_options),
        "pressure_smoother" => get(ENV, "PERF_PRESSURE_SMOOTHER", "spai0"),
        "system_smoother" => get(ENV, "PERF_SYSTEM_SMOOTHER", "default"),
        "linear_solver" => get(ENV, "PERF_LINEAR_SOLVER", "gmres"),
        "mode" => get(ENV, "PERF_MODE", "ka"),
        "single_precision_linear_solver" => get(ENV, "PERF_SINGLE_PRECISION_LINEAR", "false"),
        "index_type" => get(ENV, "PERF_INDEX_TYPE", "Int"),
        "relaxation" => "false",
        "reduce_memory" => get(ENV, "PERF_REDUCE_MEMORY", "true"),
        "gs_levels_prototype" => get(ENV, "PERF_GS_LEVELS", "false"),
        "pinned_cpus" => get(ENV, "PERF_PIN_CPUS", ""),
        "gs_global_chunks_prototype" => get(ENV, "PERF_GS_GLOBAL_CHUNKS", "false"),
        "julia" => string(VERSION), "label" => label,
        "package_versions" => Dict("Jutul" => string(pkgversion(Jutul)),
            "JutulDarcy" => string(pkgversion(JutulDarcy)),
            "Krylov" => string(pkgversion(Jutul.Krylov)),
            "KernelAbstractions" => string(pkgversion(Jutul.KernelExecution.KernelAbstractions)),
            "Polyester" => string(pkgversion(Jutul.Polyester))),
        "jutul_commit" => source_state["Jutul"]["commit"],
        "darcy_commit" => source_state["JutulDarcy"]["commit"],
        "source_state" => source_state,
        "samples" => samples,
        "median_seconds" => median(s["seconds"] for s in samples))
    open(prefix * ".toml", "w") do io
        TOML.print(io, result)
    end
    if get(ENV, "PERF_PROFILE", "false") == "true"
        sim, config = new_simulator(case; solver_options)
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
