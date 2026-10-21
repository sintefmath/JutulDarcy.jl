# Fixed-input operation scaling; no timestep selection or simulation controls.
using Jutul, JutulDarcy, LinearAlgebra, Serialization, Statistics, TOML, StaticArrays
include("pin_threads.jl")
const KP = Jutul.KAPreconditioners
BLAS.set_num_threads(1)

function term_measure(f::F; inner = 10) where F
    f()
    samples = [@timed(for _ in 1:inner; f(); end) for _ in 1:9]
    return Dict("seconds" => median(s.time for s in samples)/inner,
        "bytes" => minimum(s.bytes for s in samples)/inner,
        "compile_seconds" => sum(s.compile_time for s in samples))
end

function term_row(label, f; kw...)
    result = term_measure(f; kw...)
    println(label, ": ", result)
    flush(stdout)
    return merge(Dict("term" => label), result)
end

function probe_terms()
    fixture = deserialize(ARGS[1])
    cpr = fixture.preconditioner
    A = JutulDarcy.reservoir_jacobian(fixture.linear_system)
    scalar_rhs = copy(cpr.storage.r_ps)
    scalar_rhs ./= norm(scalar_rhs)
    block_size = size(eltype(A), 1)
    rhs = collect(reinterpret(SVector{block_size, Float64}, scalar_rhs))
    x = similar(rhs)
    scalar_x = similar(scalar_rhs)
    pA = cpr.pressure_precond.factor.levels[1].A
    prhs = copy(cpr.storage.r_p)
    prhs ./= norm(prhs)
    px = similar(prhs)
    context = KernelAbstractionsContext(Jutul.KernelExecution.KernelAbstractions.CPU())
    # The fixture's CSR worker metadata comes from its creating process. Use the
    # same arrays with current-process metadata for this native CSR product.
    A = Jutul.StaticSparsityMatrixCSR(A.nzval, A.colval, A.rowptr,
        size(A, 1), size(A, 2), A.backend;
        nthreads = Threads.nthreads(:default), minbatch = A.minbatch, thread_type = :batch)
    results = Any[]
    push!(results, term_row("block_csr_product", () -> mul!(x, A, rhs)))
    # Keep both dispatches explicit. The host three-argument product and AMG
    # residual use runtime batches even when the CSR's stored worker count is 1.
    push!(results, term_row("native_pressure_product_metadata_$(pA.nthreads)", () -> mul!(px, pA, prhs)))
    push!(results, term_row("amg_pressure_residual", () -> KP.residual!(px, pA, prhs, prhs)))
    push!(results, term_row("system_apply", () -> Jutul.apply!(scalar_x, cpr.system_precond, scalar_rhs)))
    push!(results, term_row("system_numeric_update", () ->
        Jutul.update_preconditioner!(cpr.system_precond, A, scalar_rhs, context, nothing); inner = 3))
    push!(results, term_row("cpr_apply_spai_pressure", () -> Jutul.apply!(scalar_x, cpr, scalar_rhs)))
    push!(results, term_row("krylov_dot", () -> dot(scalar_rhs, scalar_rhs); inner = 50))
    push!(results, term_row("krylov_norm", () -> norm(scalar_rhs); inner = 50))
    push!(results, term_row("krylov_axpy", () -> axpy!(0.01, scalar_rhs, scalar_x); inner = 50))
    pressure_results = Any[]
    for (method, steps) in ((:spai0, 1), (:spai0, 2), (:gauss_seidel, 1))
        label = "$(method)_steps$(steps)"
        pressure = AMGPreconditioner(:hmis; smoother_type = method, npre = steps)
        Jutul.update_preconditioner!(pressure, pA, prhs, context, nothing)
        H = pressure.factor
        row = Dict{String, Any}("smoother" => label, "terms" => Any[], "levels" => Any[])
        push!(row["terms"], term_row("pressure_apply_$(label)", () -> Jutul.apply!(px, pressure, prhs)))
        push!(row["terms"], term_row("pressure_numeric_update_$(label)", () ->
            Jutul.partial_update_preconditioner!(pressure, pA, prhs, context, nothing); inner = 3))
        Jutul.apply!(px, pressure, prhs)
        for (index, level) in enumerate(H.levels)
            index == length(H.levels) && break
            b = index == 1 ? prhs : H.levels[index - 1].rhs
            iterate = zeros(eltype(b), size(level.A, 1))
            operations = Any[]
            push!(operations, term_row("level$(index)_down_$(label)", () ->
                KP.smooth_level!(iterate, level.A, b, level.smoother, H.options.smoother.steps;
                    residual = b, zero_initial = true)))
            push!(operations, term_row("level$(index)_up_$(label)", () ->
                KP.smooth_result!(iterate, level.A, b, level.smoother, H.options.smoother.steps)))
            push!(operations, term_row("level$(index)_residual_$(label)", () ->
                KP.residual!(level.residual, level.A, iterate, b)))
            push!(operations, term_row("level$(index)_restrict_$(label)", () ->
                KP.restrict!(level.rhs, level.Pt, level.P, level.residual, H.backend, H.block_size)))
            push!(operations, term_row("level$(index)_prolong_$(label)", () ->
                KP.prolong!(iterate, level.P, level.correction, H.backend, H.block_size)))
            push!(row["levels"], Dict("index" => index, "rows" => size(level.A, 1),
                "nonzeros" => length(level.A.nzval), "terms" => operations))
        end
        push!(pressure_results, row)
    end
    result = Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
        "julia" => string(VERSION), "blas_threads" => BLAS.get_num_threads(),
        "pinned_cpus" => get(ENV, "PERF_PIN_CPUS", ""),
        "cells" => size(A, 1), "block_size" => block_size,
        "terms" => results, "pressure" => pressure_results,
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`))
    open(ARGS[2], "w") do io
        TOML.print(io, result)
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    probe_terms()
end
