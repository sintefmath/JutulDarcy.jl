# Use the existing native triangular row kernel and split CSR construction with
# KA's dependency schedule. Repack saved factors; do not refactor another iterate.
include("probe_terms.jl")
using Polyester, Random
const SC = Jutul.StaticCSR

function split_ilu_levels!(x, work, L, U, diagonal, state, b, batch_size)
    lower_rows, lower_offsets = state.factor_rows, state.factor_offsets
    for level in 1:length(lower_offsets)-1
        lower_start = lower_offsets[level]
        count = lower_offsets[level + 1] - lower_start
        @batch minbatch=batch_size for q in 1:count
            i = lower_rows[lower_start + q - 1]
            SC.invert_row!(work, L, nothing, b, i, i)
        end
    end
    upper_rows, upper_offsets = state.upper_rows, state.upper_offsets
    for level in 1:length(upper_offsets)-1
        upper_start = upper_offsets[level]
        count = upper_offsets[level + 1] - upper_start
        @batch minbatch=batch_size for q in 1:count
            i = upper_rows[upper_start + q - 1]
            SC.invert_row!(x, U, diagonal, work, i, i)
        end
    end
    return x
end

function probe_ilu_storage()
    state = deserialize(ARGS[1]).preconditioner.system_precond.factor
    n = state.n
    @assert isone(KP.smoother_damping(state))
    matrix = SC.StaticSparsityMatrixCSR(state.factors, state.colval, state.rowptr,
        n, n, state.backend; nthreads = 1, minbatch = state.block_size)
    L, lower_map = SC.fixed_block(matrix; lower = true)
    U, upper_map = SC.fixed_block(matrix; lower = false)
    T = eltype(state.work)
    Random.seed!(1234)
    rhs, x = randn(T, n), similar(state.work)
    KP.ilu_solve_cpu!(x, state, rhs, state.work)
    reference = copy(x)
    rows = Any[term_row("original", () -> KP.ilu_solve_cpu!(x, state, rhs, state.work))]
    for batch_size in (1, 8, 32, 128)
        f = () -> split_ilu_levels!(x, state.work, L, U, state.inverse_diagonal, state, rhs, batch_size)
        f()
        scalar_type = KP.matrix_scalar_type(T)
        relative_error = norm(reinterpret(scalar_type, x) - reinterpret(scalar_type, reference))/norm(reinterpret(scalar_type, reference))
        @assert relative_error < 1e-12
        push!(rows, merge(Dict("relative_error" => relative_error,
            "bitwise_equal" => isequal(x, reference)), term_row("split_minbatch_$batch_size", f)))
    end
    serial() = SC.ilu_solve!(x, L, U, state.inverse_diagonal, rhs)
    push!(rows, term_row("split_native_serial", serial))
    push!(rows, term_row("refresh_split_values", () -> begin
        SC.inner_update!(L.nzval, state.factors, lower_map)
        SC.inner_update!(U.nzval, state.factors, upper_map)
    end))
    push!(rows, term_row("original_end", () -> KP.ilu_solve_cpu!(x, state, rhs, state.work)))
    result = Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
        "rows" => n, "terms" => rows, "source_nonzeros" => length(state.factors),
        "split_nonzeros" => length(L.nzval) + length(U.nzval),
        "extra_storage_bytes" => Base.summarysize((L, U, lower_map, upper_map)),
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`))
    open(ARGS[2], "w") do io
        TOML.print(io, result)
    end
end
probe_ilu_storage()
