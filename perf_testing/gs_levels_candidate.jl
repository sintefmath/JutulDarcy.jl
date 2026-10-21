# Process-local experiment, included into KAPreconditioners by the replay only.
# Production would keep this cache on the smoother state.
include(joinpath(@__DIR__, "level_rows_candidate.jl"))
struct PerfGSLevels{Ti}
    rowptr::Vector{Ti}
    colval::Vector{Ti}
    offsets::Vector{Int}
    rows::Vector{Int}
    width::Int
    source_rowptr::Vector{Ti}
    source_colval::Vector{Ti}
end
const perf_gs_cache = Dict{UInt, Tuple{WeakRef, Any}}()
const perf_gs_use_levels = Ref(false)
const perf_gs_minimum_rows = Ref(16_384)
const perf_gs_minbatch = Ref(8)
const perf_gs_global_chunks = Ref(false)
const perf_gs_raw = Ref(parse(Bool, get(ENV, "PERF_GS_RAW", "false")))

@inline function perf_gs_raw_row!(x, target, rowptr, colval, matrix_values, b,
        inverse_diagonal, damping, i, ::Val{correction}) where correction
    @inbounds begin
        value = b[i]
        @simd for k in rowptr[i]:rowptr[i + 1]-1
            value -= matrix_values[k]*target[colval[k]]
        end
        delta = damping*inverse_diagonal[i]*value
        if correction
            target[i] = delta
        end
        x[i] += delta
    end
    return nothing
end

function perf_gs_raw_levels!(x, target, A, b, inverse_diagonal, damping,
        schedule::PerfGSLevels, forward, ::Val{correction}) where correction
    correction && fill!(target, zero(eltype(target)))
    rowptr, colval, matrix_values = A.rowptr, A.colval, A.nzval
    offsets, scheduled_rows = schedule.offsets, schedule.rows
    nlevels = length(offsets) - 1
    batch_size = perf_gs_minbatch[]
    workers, fixed = Threads.nthreads(:default), perf_gs_global_chunks[]
    mode = Val(correction)
    for batch in 1:nlevels
        level = forward ? batch : nlevels - batch + 1
        level_start = offsets[level]
        count = offsets[level + 1] - level_start
        if fixed && count >= 2*workers
            @batch for worker in 1:workers
                low, high = perf_level_interval(scheduled_rows, level_start, count, worker, length(x), workers)
                @inbounds for q in low:high
                    perf_gs_raw_row!(x, target, rowptr, colval, matrix_values, b,
                        inverse_diagonal, damping, scheduled_rows[q], mode)
                end
            end
        elseif fixed
            @inbounds for q in level_start:level_start+count-1
                perf_gs_raw_row!(x, target, rowptr, colval, matrix_values, b,
                    inverse_diagonal, damping, scheduled_rows[q], mode)
            end
        else
            @batch minbatch=batch_size for q in 1:count
                @inbounds i = scheduled_rows[level_start + q - 1]
                perf_gs_raw_row!(x, target, rowptr, colval, matrix_values, b,
                    inverse_diagonal, damping, i, mode)
            end
        end
    end
    return x
end

function perf_gs_schedule(A, inverse_diagonal; validate = true)
    n, nnz = matrix_nrows(A), matrix_nonzeros(A)
    entry = get(perf_gs_cache, objectid(inverse_diagonal), nothing)
    saved = isnothing(entry) || entry[1].value !== inverse_diagonal ? nothing : entry[2]
    if !isnothing(saved) && length(saved.rowptr) == n + 1 && length(saved.colval) == nnz &&
            saved.source_rowptr === A.rowptr && saved.source_colval === A.colval &&
            (!validate || (same_prefix(A.rowptr, saved.rowptr, n + 1) && same_prefix(A.colval, saved.colval, nnz)))
        return saved
    end
    levels = ones(Int, n)
    @inbounds for i in 1:n
        level = levels[i]
        for k in A.rowptr[i]:(A.rowptr[i + 1] - 1)
            j = A.colval[k]
            j < i && (level = max(level, levels[j] + 1))
        end
        levels[i] = level
        # Include upper-only edges as ordering constraints as well. They must
        # still read the old higher-index value during a forward sweep.
        for k in A.rowptr[i]:(A.rowptr[i + 1] - 1)
            j = A.colval[k]
            j > i && (levels[j] = max(levels[j], level + 1))
        end
    end
    counts = zeros(Int, maximum(levels; init = 0))
    for level in levels
        counts[level] += 1
    end
    offsets = cumsum(vcat(1, counts))
    cursor = copy(offsets)
    rows = Vector{Int}(undef, n)
    for i in 1:n
        level = levels[i]
        rows[cursor[level]] = i
        cursor[level] += 1
    end
    saved = PerfGSLevels(copy(view(A.rowptr, 1:(n + 1))), copy(view(A.colval, 1:nnz)),
        offsets, rows, maximum(counts; init = 0), A.rowptr, A.colval)
    for (key, (reference, _)) in perf_gs_cache
        isnothing(reference.value) && delete!(perf_gs_cache, key)
    end
    perf_gs_cache[objectid(inverse_diagonal)] = (WeakRef(inverse_diagonal), saved)
    return saved
end

function update_smoother!(state::GaussSeidelState{D},
        A::StaticSparsityMatrixCSR{Tv, Ti, <:Vector, <:Vector, <:Vector}) where {D <: Vector, Tv, Ti}
    invoke(update_smoother!, Tuple{GaussSeidelState, typeof(A)}, state, A)
    if perf_gs_use_levels[] && Threads.nthreads(:default) > 1 && matrix_nrows(A) >= perf_gs_minimum_rows[]
        perf_gs_schedule(A, state.inverse_diagonal)
    end
    return state
end

function perf_gs_enabled(A, rows)
    n = matrix_nrows(A)
    return perf_gs_use_levels[] && Threads.nthreads(:default) > 1 &&
        matrix_backend(A) isa KernelAbstractions.CPU && n >= perf_gs_minimum_rows[] &&
        length(rows) == n && abs(step(rows)) == 1 &&
        ((first(rows) == 1 && last(rows) == n) || (first(rows) == n && last(rows) == 1))
end

function perf_gs_sweep_levels!(x, A, b, inverse_diagonal, damping, schedule::PerfGSLevels, forward)
    if perf_gs_raw[]
        return perf_gs_raw_levels!(x, x, A, b, inverse_diagonal, damping, schedule, forward, Val(false))
    end
    offsets, scheduled_rows = schedule.offsets, schedule.rows
    nlevels = length(offsets) - 1
    operation = row -> gauss_seidel_sweep!(x, A, b, inverse_diagonal, damping, (row,))
    for batch in 1:nlevels
        level = forward ? batch : nlevels - batch + 1
        start = offsets[level]
        count = offsets[level + 1] - start
        perf_level_rows!(operation, scheduled_rows, start, count, length(x),
            perf_gs_minbatch[], perf_gs_global_chunks[])
    end
    return x
end

function gauss_seidel_sweep!(x::Vector, A::StaticSparsityMatrixCSR{Tv, Ti, <:Vector, <:Vector, <:Vector},
        b::Vector, inverse_diagonal::Vector, damping, rows::AbstractRange{<:Integer}) where {Tv, Ti}
    if perf_gs_enabled(A, rows)
        schedule = perf_gs_schedule(A, inverse_diagonal; validate = false)
        if schedule.width > 1
            return perf_gs_sweep_levels!(x, A, b, inverse_diagonal, damping, schedule, step(rows) > 0)
        end
    end
    return invoke(gauss_seidel_sweep!, Tuple{Any, Any, Any, Any, Any, Any},
        x, A, b, inverse_diagonal, damping, rows)
end

function perf_gs_correction_levels!(x, correction, A, residual, inverse_diagonal, damping,
        schedule::PerfGSLevels, forward)
    if perf_gs_raw[]
        return perf_gs_raw_levels!(x, correction, A, residual, inverse_diagonal, damping, schedule, forward, Val(true))
    end
    fill!(correction, zero(eltype(correction)))
    offsets, scheduled_rows = schedule.offsets, schedule.rows
    nlevels = length(offsets) - 1
    operation = function(i)
        @inbounds begin
            value = residual[i]
            @simd for k in A.rowptr[i]:(A.rowptr[i + 1] - 1)
                value -= A.nzval[k] * correction[A.colval[k]]
            end
            correction[i] = damping * inverse_diagonal[i] * value
            x[i] += correction[i]
        end
    end
    for batch in 1:nlevels
        level = forward ? batch : nlevels - batch + 1
        start = offsets[level]
        count = offsets[level + 1] - start
        perf_level_rows!(operation, scheduled_rows, start, count, length(x),
            perf_gs_minbatch[], perf_gs_global_chunks[])
    end
    return x
end

function gauss_seidel_correction!(x::Vector, correction::Vector,
        A::StaticSparsityMatrixCSR{Tv, Ti, <:Vector, <:Vector, <:Vector}, residual::Vector,
        inverse_diagonal::Vector, damping, rows::AbstractRange{<:Integer}) where {Tv, Ti}
    if perf_gs_enabled(A, rows)
        schedule = perf_gs_schedule(A, inverse_diagonal; validate = false)
        if schedule.width > 1
            return perf_gs_correction_levels!(x, correction, A, residual,
                inverse_diagonal, damping, schedule, step(rows) > 0)
        end
    end
    return invoke(gauss_seidel_correction!, Tuple{Any, Any, Any, Any, Any, Any, Any},
        x, correction, A, residual, inverse_diagonal, damping, rows)
end
