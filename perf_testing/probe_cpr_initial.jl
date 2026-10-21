# Test whether the first system solve can write directly to the CPR output.
using Jutul, JutulDarcy, LinearAlgebra, Serialization, Statistics, TOML

function direct_initial_apply!(x, cpr, rhs)
    storage = cpr.storage
    r, buffer, A = storage.r_ps, storage.x_ps, storage.A_ps
    copyto!(r, rhs)
    if cpr.npre > 0
        Jutul.apply!(x, cpr.system_precond, r)
        JutulDarcy.correct_residual!(r, A, x)
        if cpr.npre > 1
            JutulDarcy.apply_cpr_smoother!(x, r, buffer, cpr.system_precond, A, cpr.npre - 1)
        end
    else
        fill!(x, zero(eltype(x)))
    end
    JutulDarcy.apply_cpr_pressure_stage!(cpr, storage, r)
    if cpr.npost > 0
        JutulDarcy.correct_residual_and_increment_pressure!(r, x, storage.p,
            storage.block_size, buffer, A, storage.p_buffer)
        JutulDarcy.apply_cpr_smoother!(x, r, buffer, cpr.system_precond, A, cpr.npost; skip_last = true)
    else
        JutulDarcy.increment_pressure!(x, storage.p, storage.block_size, storage.p_buffer)
    end
    return x
end

function application_trial(f)
    f()
    samples = [@timed(f()) for _ in 1:9]
    return Dict("seconds" => median(s.time for s in samples),
        "bytes" => minimum(s.bytes for s in samples))
end

function probe_cpr_initial()
    fixture = deserialize(ARGS[1])
    cpr = fixture.preconditioner
    matrix = JutulDarcy.reservoir_jacobian(fixture.linear_system)
    symbolic = Jutul.KAPreconditioners.ilu_symbolic(matrix)
    symbolic_setup = application_trial(() -> Jutul.KAPreconditioners.ilu_symbolic(matrix))
    transpose_setup = application_trial(() -> Jutul.KAPreconditioners.transpose_positions(
        symbolic.rowptr, symbolic.colval, size(matrix, 1)))
    println("ILU symbolic setup: ", (symbolic_setup, transpose_setup))
    flush(stdout)
    rhs = copy(cpr.storage.r_ps)
    rhs ./= norm(rhs)
    x = similar(rhs)
    fixture = nothing
    GC.gc()
    results = Any[]
    for pre in 0:2, post in 0:1
        cpr.npre, cpr.npost = pre, post
        fill!(x, NaN)
        Jutul.apply!(x, cpr, rhs)
        reference = copy(x)
        fill!(x, NaN)
        direct_initial_apply!(x, cpr, rhs)
        @assert x ≈ reference
        row = Dict("npre" => pre, "npost" => post, "bitwise_equal" => x == reference,
            "original" => application_trial(() -> Jutul.apply!(x, cpr, rhs)),
            "direct_initial" => application_trial(() -> direct_initial_apply!(x, cpr, rhs)))
        push!(results, row)
        println(row)
        flush(stdout)
    end
    open(ARGS[2], "w") do io
        TOML.print(io, Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
            "symbolic_setup" => symbolic_setup, "transpose_setup" => transpose_setup,
            "results" => results))
    end
end
probe_cpr_initial()
