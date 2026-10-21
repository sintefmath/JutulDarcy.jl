# Compare triangular traversal on a seven-point block-CSR matrix using the
# existing level schedule and stored diagonal. No runtime method replacement.
using Jutul, Polyester, StaticArrays, SparseArrays, LinearAlgebra, Statistics, TOML, Random
const KP = Jutul.KAPreconditioners
Random.seed!(1234)

function triangular_ranges!(x, state, b)
    factors, diagonal = state.factors, state.diagonal_positions
    rowptr, colval = state.rowptr, state.colval
    work = state.work
    lower_rows = state.factor_rows
    lower_offsets = state.factor_offsets
    for level in 1:(length(lower_offsets) - 1)
        lower_start = lower_offsets[level]
        count = lower_offsets[level + 1] - lower_start
        @batch for q in 1:count
            @inbounds begin
                i = lower_rows[lower_start + q - 1]
                value = b[i]
                for k in rowptr[i]:(diagonal[i] - 1)
                    value -= factors[k] * work[colval[k]]
                end
                work[i] = value
            end
        end
    end
    upper_rows = state.upper_rows
    upper_offsets = state.upper_offsets
    damping = KP.smoother_damping(state)
    inverse_diagonal = state.inverse_diagonal
    for level in 1:(length(upper_offsets) - 1)
        upper_start = upper_offsets[level]
        count = upper_offsets[level + 1] - upper_start
        @batch for q in 1:count
            @inbounds begin
                i = upper_rows[upper_start + q - 1]
                value = work[i]
                for k in (diagonal[i] + 1):(rowptr[i + 1] - 1)
                    value -= factors[k] * x[colval[k]]
                end
                x[i] = damping * (inverse_diagonal[i] * value)
            end
        end
    end
    return x
end

function trial(f; repeats=10)
    f()
    samples = [@timed(for _ in 1:repeats; f(); end) for _ in 1:5]
    return Dict("seconds" => median(s.time for s in samples)/repeats,
        "bytes" => median(s.bytes for s in samples)/repeats)
end

function main()
    n = 192000
    scalar = spdiagm(-4800 => fill(-0.8, n-4800), -80 => fill(-0.9, n-80),
        -1 => fill(-1.1, n-1), 0 => fill(8.0, n),
        1 => fill(-1.0, n-1), 80 => fill(-0.85, n-80), 4800 => fill(-0.7, n-4800))
    # Remove wraparound edges at Cartesian row/plane boundaries. Retaining them
    # would turn the lower graph into a chain with one row per dependency level.
    for column in 1:n, position in nzrange(scalar, column)
        row = scalar.rowval[position]
        offset = abs(row-column)
        if (offset == 1 && (row-1) ÷ 80 != (column-1) ÷ 80) ||
                (offset == 80 && (row-1) ÷ 4800 != (column-1) ÷ 4800)
            scalar.nzval[position] = 0
        end
    end
    dropzeros!(scalar)
    results = []
    for block in (1, 2, 3)
        if block == 1
            matrix = KP.csr_matrix(scalar)
            rhs = randn(n)
        else
            value_type = SMatrix{block, block, Float64, block*block}
            blocks = map(v -> value_type(v*I), scalar.nzval)
            matrix = KP.csr_matrix(SparseMatrixCSC(n, n, scalar.colptr, scalar.rowval, blocks))
            rhs = randn(SVector{block,Float64}, n)
        end
        state = Jutul.setup_ka_smoother(matrix, KP.ILU0())
        x = similar(rhs)
        original = () -> KP.ilu_solve_cpu!(x, state, rhs, state.work)
        ranged = () -> triangular_ranges!(x, state, rhs)
        reference = copy(original())
        @assert ranged() == reference
        row = Dict("block_size" => block, "rows" => n,
            "threads" => Threads.nthreads(), "levels" => length(state.factor_offsets)-1,
            "original" => trial(original), "ranged" => trial(ranged))
        push!(results, row)
        println(row)
        flush(stdout)
    end
    label = get(ENV, "PERF_LABEL", "ilu_grid_ranges")
    open(joinpath(@__DIR__, "results", "$(label)_t$(Threads.nthreads()).toml"), "w") do io
        TOML.print(io, Dict("results" => results))
    end
end
main()
