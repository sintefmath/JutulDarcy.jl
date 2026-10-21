# Separate dependency-level dispatch cost from traversal locality. Same factors,
# dependency schedule and within-row arithmetic; no simulation controls change.
include("probe_terms.jl")
using Polyester, Random

function ilu_minbatch!(x, state, b, batch_size)
    factors, diagonal = state.factors, state.inverse_diagonal
    rowptr, colval, work = state.rowptr, state.colval, state.work
    lower_rows, lower_offsets = state.factor_rows, state.factor_offsets
    for level in 1:length(lower_offsets)-1
        lower_start = lower_offsets[level]
        count = lower_offsets[level + 1] - lower_start
        @batch minbatch=batch_size for q in 1:count
            @inbounds begin
                i = lower_rows[lower_start + q - 1]
                value = b[i]
                for k in rowptr[i]:rowptr[i + 1]-1
                    j = colval[k]
                    j < i && (value -= factors[k]*work[j])
                end
                work[i] = value
            end
        end
    end
    upper_rows, upper_offsets = state.upper_rows, state.upper_offsets
    damping = KP.smoother_damping(state)
    for level in 1:length(upper_offsets)-1
        upper_start = upper_offsets[level]
        count = upper_offsets[level + 1] - upper_start
        @batch minbatch=batch_size for q in 1:count
            @inbounds begin
                i = upper_rows[upper_start + q - 1]
                value = work[i]
                for k in rowptr[i]:rowptr[i + 1]-1
                    j = colval[k]
                    j > i && (value -= factors[k]*x[j])
                end
                x[i] = damping*(diagonal[i]*value)
            end
        end
    end
    return x
end

function ilu_natural!(x, state, b)
    factors, diagonal = state.factors, state.inverse_diagonal
    rowptr, colval, work = state.rowptr, state.colval, state.work
    @inbounds for i in eachindex(b)
        value = b[i]
        for k in rowptr[i]:rowptr[i + 1]-1
            j = colval[k]
            j < i && (value -= factors[k]*work[j])
        end
        work[i] = value
    end
    damping = KP.smoother_damping(state)
    @inbounds for i in reverse(eachindex(b))
        value = work[i]
        for k in rowptr[i]:rowptr[i + 1]-1
            j = colval[k]
            j > i && (value -= factors[k]*x[j])
        end
        x[i] = damping*(diagonal[i]*value)
    end
    return x
end

function level_store!(scratch, offsets, batch_size)
    for level in 1:length(offsets)-1
        count = offsets[level + 1] - offsets[level]
        @batch minbatch=batch_size for q in 1:count
            @inbounds scratch[q] = q
        end
    end
    return nothing
end

function probe_ilu_dispatch()
    state = deserialize(ARGS[1]).preconditioner.system_precond.factor
    T = eltype(state.work)
    Random.seed!(1234)
    rhs = randn(T, length(state.work))
    x = similar(rhs)
    KP.ilu_solve_cpu!(x, state, rhs, state.work)
    reference = copy(x)
    rows = Any[term_row("original", () -> KP.ilu_solve_cpu!(x, state, rhs, state.work))]
    for size in (1, 8, 32, 128, 512, 2048, 8192)
        f = () -> ilu_minbatch!(x, state, rhs, size)
        f()
        @assert isequal(x, reference)
        push!(rows, merge(Dict("bitwise_equal" => true), term_row("minbatch_$size", f)))
    end
    ilu_natural!(x, state, rhs)
    scalar_type = KP.matrix_scalar_type(T)
    scalar_x, scalar_ref = reinterpret(scalar_type, x), reinterpret(scalar_type, reference)
    relative_error = norm(scalar_x - scalar_ref)/norm(scalar_ref)
    @assert relative_error < 1e-12
    push!(rows, merge(Dict("bitwise_equal" => isequal(x, reference),
        "relative_error" => relative_error),
        term_row("natural_serial", () -> ilu_natural!(x, state, rhs))))
    counts = vcat(diff(state.factor_offsets), diff(state.upper_offsets))
    scratch = zeros(Int, maximum(counts))
    for size in (1, 32, 128, 512)
        f = () -> begin
            level_store!(scratch, state.factor_offsets, size)
            level_store!(scratch, state.upper_offsets, size)
        end
        push!(rows, term_row("level_store_minbatch_$size", f))
    end
    result = Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
        "rows" => length(rhs), "block_type" => string(T), "terms" => rows,
        "level_row_counts" => counts, "levels" => length(counts),
        "minimum_level_rows" => minimum(counts), "median_level_rows" => median(counts),
        "maximum_level_rows" => maximum(counts),
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`))
    open(ARGS[2], "w") do io
        TOML.print(io, result)
    end
end
probe_ilu_dispatch()
