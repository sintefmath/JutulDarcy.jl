# Replay the actual pressure matrix/RHS saved by PERF_SAVE_FIXTURE=true.
using Jutul, JutulDarcy, LinearAlgebra, Serialization, Statistics, TOML
include("pin_threads.jl")
Base.include(Jutul.KAPreconditioners, joinpath(@__DIR__, "gs_levels_candidate.jl"))
Base.include(Jutul.KAPreconditioners, joinpath(@__DIR__, "hybrid_smoother_candidate.jl"))

function pressure_measure(f; repeats = 9)
    f()
    samples = [@timed(f()) for _ in 1:repeats]
    return Dict("seconds" => median(s.time for s in samples),
        "bytes" => minimum(s.bytes for s in samples))
end

function probe_pressure()
    fixture = deserialize(ARGS[1])
    A = fixture.preconditioner.pressure_precond.factor.levels[1].A
    rhs = copy(fixture.preconditioner.storage.r_p)
    rhs ./= norm(rhs)
    fixture = nothing
    GC.gc()
    choices = (
        ("gs", :hmis, (smoother_type = :gauss_seidel,)),
        ("gs_levels", :hmis, (smoother_type = :gauss_seidel,)),
        ("gs_chunks", :hmis, (smoother_type = :gauss_seidel,)),
        ("hybrid8k", :hmis, (smoother = Jutul.KAPreconditioners.PerfHybridSmoother(8192),)),
        ("hybrid32k", :hmis, (smoother = Jutul.KAPreconditioners.PerfHybridSmoother(32768),)),
        ("hybrid32k_2", :hmis, (smoother = Jutul.KAPreconditioners.PerfHybridSmoother(32768; steps = 2),)),
        ("hybrid128k", :hmis, (smoother = Jutul.KAPreconditioners.PerfHybridSmoother(131072),)),
        ("hybrid256k", :hmis, (smoother = Jutul.KAPreconditioners.PerfHybridSmoother(262144),)),
        ("hybrid256k_2", :hmis, (smoother = Jutul.KAPreconditioners.PerfHybridSmoother(262144; steps = 2),)),
        ("spai1", :hmis, (smoother_type = :spai0,)),
        ("spai2", :hmis, (smoother_type = :spai0, npre = 2)),
        ("spai4", :hmis, (smoother_type = :spai0, npre = 4)),
        ("ilu1", :hmis, (smoother_type = :ilu0,)),
        ("ilu07", :hmis, (smoother_type = :ilu0, damping = 0.7)),
        ("dilu1", :hmis, (smoother_type = :dilu,)),
        ("dilu07", :hmis, (smoother_type = :dilu, damping = 0.7)),
        ("gs_coarse32", :hmis, (smoother_type = :gauss_seidel, coarse_size = 32)),
        ("gs_aggressive1", :hmis, (smoother_type = :gauss_seidel, aggressive_levels = 1)),
        ("gs_theta025", :hmis, (smoother_type = :gauss_seidel, theta = 0.25)),
        ("gs_theta01", :hmis, (smoother_type = :gauss_seidel, theta = 0.1)),
        ("spai2_theta025", :hmis, (smoother_type = :spai0, npre = 2, theta = 0.25)),
        ("rs_gs", :ruge_stuben, (smoother_type = :gauss_seidel,)),
        ("aggregation_gs", :aggregation, (smoother_type = :gauss_seidel,)),
    )
    selected = split(get(ENV, "PERF_PRESSURE_POLICIES", join(first.(choices), " ")))
    context = KernelAbstractionsContext(Jutul.KernelExecution.KernelAbstractions.CPU())
    output = Any[]
    for (label, method, options) in choices
        label in selected || continue
        println("Pressure replay: ", label)
        flush(stdout)
        Jutul.KAPreconditioners.perf_gs_use_levels[] = label in ("gs_levels", "gs_chunks")
        Jutul.KAPreconditioners.perf_gs_global_chunks[] = label == "gs_chunks"
        result = try
            pressure = AMGPreconditioner(method; options...)
            Jutul.update_preconditioner!(pressure, A, rhs, context, nothing)
            x = similar(rhs)
            application = pressure_measure(() -> Jutul.apply!(x, pressure, rhs))
            numeric_update = pressure_measure(() -> Jutul.partial_update_preconditioner!(pressure, A, rhs, context, nothing); repeats = 3)
            solver = Jutul.Krylov.GmresSolver(A, rhs, 50)
            solve_pressure() = Jutul.Krylov.gmres!(solver, A, rhs;
                N = pressure.factor, rtol = 0.005, atol = 0.0, itmax = 100, history = true)
            solve_pressure()
            solve = @timed solve_pressure()
            levels = Any[]
            for level in pressure.factor.levels
                smoother = level.smoother
                row = Dict{String, Any}("rows" => size(level.A, 1), "nonzeros" => length(level.A.nzval))
                if hasproperty(smoother, :inverse_diagonal)
                    row["negative_inverse_diagonals"] = count(<(0), smoother.inverse_diagonal)
                    row["nonfinite_inverse_diagonals"] = count(!isfinite, smoother.inverse_diagonal)
                    if label in ("gs_levels", "gs_chunks")
                        schedule = Jutul.KAPreconditioners.perf_gs_schedule(level.A, smoother.inverse_diagonal)
                        row["gs_levels"] = length(schedule.offsets) - 1
                        row["gs_maximum_width"] = schedule.width
                    end
                end
                push!(levels, row)
            end
            Dict("label" => label, "method" => String(method), "options" => repr(options),
                "apply" => application, "numeric_update" => numeric_update,
                "pressure_gmres_seconds" => solve.time, "iterations" => solver.stats.niter,
                "solved" => solver.stats.solved,
                "true_relative_residual" => norm(rhs - A*solver.x)/norm(rhs), "levels" => levels)
        catch err
            err isa Union{LinearAlgebra.SingularException, LinearAlgebra.PosDefException} || rethrow()
            Dict("label" => label, "method" => String(method), "options" => repr(options),
                "status" => "failed", "error" => sprint(showerror, err))
        end
        push!(output, result)
        println(result)
        flush(stdout)
        open(ARGS[2], "w") do io
            TOML.print(io, Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
                "julia" => string(VERSION), "pinned_cpus" => get(ENV, "PERF_PIN_CPUS", ""),
                "results" => output))
        end
        pressure = solver = nothing
        GC.gc()
    end
end
probe_pressure()
