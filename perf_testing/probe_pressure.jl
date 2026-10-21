# Replay the actual pressure matrix/RHS saved by PERF_SAVE_FIXTURE=true.
using Jutul, JutulDarcy, LinearAlgebra, Serialization, Statistics, TOML

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
        "gs" => (smoother_type = :gauss_seidel,),
        "spai1" => (smoother_type = :spai0,),
        "spai2" => (smoother_type = :spai0, npre = 2),
        "spai4" => (smoother_type = :spai0, npre = 4),
        "ilu1" => (smoother_type = :ilu0,),
        "ilu07" => (smoother_type = :ilu0, damping = 0.7),
        "dilu1" => (smoother_type = :dilu,),
        "dilu07" => (smoother_type = :dilu, damping = 0.7),
        "gs_coarse32" => (smoother_type = :gauss_seidel, coarse_size = 32),
        "gs_aggressive1" => (smoother_type = :gauss_seidel, aggressive_levels = 1),
    )
    selected = split(get(ENV, "PERF_PRESSURE_POLICIES", join(first.(choices), " ")))
    context = KernelAbstractionsContext(Jutul.KernelExecution.KernelAbstractions.CPU())
    output = Any[]
    for (label, options) in choices
        label in selected || continue
        println("Pressure replay: ", label)
        flush(stdout)
        pressure = AMGPreconditioner(:hmis; options...)
        Jutul.update_preconditioner!(pressure, A, rhs, context, nothing)
        x = similar(rhs)
        application = pressure_measure(() -> Jutul.apply!(x, pressure, rhs))
        numeric_update = pressure_measure(() -> Jutul.partial_update_preconditioner!(pressure, A, rhs, context, nothing); repeats = 3)
        solver = Jutul.Krylov.GmresSolver(A, rhs; memory = 50)
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
            end
            push!(levels, row)
        end
        row = Dict("label" => label, "options" => repr(options),
            "apply" => application, "numeric_update" => numeric_update,
            "pressure_gmres_seconds" => solve.time, "iterations" => solver.stats.niter,
            "solved" => solver.stats.solved,
            "true_relative_residual" => norm(rhs - A*solver.x)/norm(rhs), "levels" => levels)
        push!(output, row)
        println(row)
        flush(stdout)
        open(ARGS[2], "w") do io
            TOML.print(io, Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
                "julia" => string(VERSION), "results" => output))
        end
        pressure = solver = nothing
        GC.gc()
    end
end
probe_pressure()
