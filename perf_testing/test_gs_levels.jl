using Jutul, Test, SparseArrays, LinearAlgebra, Random
const KP = Jutul.KAPreconditioners
Base.include(KP, joinpath(@__DIR__, "gs_levels_candidate.jl"))
KP.perf_gs_use_levels[] = true
KP.perf_gs_minimum_rows[] = 0
KP.perf_gs_minbatch[] = 1

function check_schedule(A, schedule)
    levels = zeros(Int, size(A, 1))
    for level in 1:(length(schedule.offsets)-1)
        for q in schedule.offsets[level]:(schedule.offsets[level+1]-1)
            levels[schedule.rows[q]] = level
        end
    end
    for i in axes(A, 1), k in A.rowptr[i]:(A.rowptr[i+1]-1)
        j = A.colval[k]
        if i < j
            @test levels[i] < levels[j]
        elseif i > j
            @test levels[j] < levels[i]
        end
    end
end

@testset "Natural-order GS dependency schedule" begin
    rng = MersenneTwister(7)
    for fixed in (false, true), T in (Float32, Float64), pattern in 1:12
        KP.perf_gs_global_chunks[] = fixed
        M = SparseMatrixCSC{T, Int}(sprand(rng, 16, 16, 0.15))
        if pattern % 3 == 0
            M = triu(M, 1)
        elseif pattern % 3 == 1
            M = tril(M, -1)
        end
        M = -M + T(5)*I
        A = KP.csr_matrix(M)
        rhs, initial = randn(rng, T, 16), randn(rng, T, 16)
        state = KP.setup_smoother(A, KP.GaussSeidel())
        schedule = KP.perf_gs_schedule(A, state.inverse_diagonal)
        check_schedule(A, schedule)
        for forward in (true, false), damping in (T(0.4), one(T))
            rows = forward ? (1:16) : (16:-1:1)
            reference = copy(initial)
            invoke(KP.gauss_seidel_sweep!, Tuple{Any, Any, Any, Any, Any, Any},
                reference, A, rhs, state.inverse_diagonal, damping, rows)
            candidate = copy(initial)
            KP.perf_gs_sweep_levels!(candidate, A, rhs, state.inverse_diagonal, damping, schedule, forward)
            @test candidate == reference
            candidate = copy(initial)
            KP.gauss_seidel_sweep!(candidate, A, rhs, state.inverse_diagonal, damping, rows)
            @test candidate == reference

            reference = copy(initial)
            invoke(KP.gauss_seidel_correction!, Tuple{Any, Any, Any, Any, Any, Any, Any},
                reference, state.correction, A, rhs, state.inverse_diagonal, damping, rows)
            candidate = copy(initial)
            KP.perf_gs_correction_levels!(candidate, state.correction, A, rhs,
                state.inverse_diagonal, damping, schedule, forward)
            @test candidate == reference
            candidate = copy(initial)
            KP.gauss_seidel_correction!(candidate, state.correction, A, rhs, state.inverse_diagonal, damping, rows)
            @test candidate == reference
        end
        state.inverse_diagonal .*= 2
        @test KP.perf_gs_schedule(A, state.inverse_diagonal) === schedule
        B = KP.csr_matrix(spdiagm(0 => fill(5.0, 16)))
        replacement = KP.perf_gs_schedule(B, state.inverse_diagonal)
        @test replacement !== schedule
        check_schedule(B, replacement)
    end
end

@testset "Wide GS levels with fixed row ownership" begin
    # Independent directed chains give wide levels even at 15 workers. The
    # smaller random cases above exercise the serial fallback for narrow levels.
    rng = MersenneTwister(19)
    for T in (Float32, Float64), upper in (false, true)
        n, width = 256, 64
        M = spdiagm(0 => fill(T(5), n),
            (upper ? width : -width) => -rand(rng, T, n - width))
        A = KP.csr_matrix(M)
        rhs, initial = randn(rng, T, n), randn(rng, T, n)
        state = KP.setup_smoother(A, KP.GaussSeidel())
        schedule = KP.perf_gs_schedule(A, state.inverse_diagonal)
        @test schedule.width >= width
        check_schedule(A, schedule)
        for fixed in (false, true), forward in (false, true)
            KP.perf_gs_global_chunks[] = fixed
            rows = forward ? (1:n) : (n:-1:1)
            reference, candidate = copy(initial), copy(initial)
            invoke(KP.gauss_seidel_sweep!, Tuple{Any, Any, Any, Any, Any, Any},
                reference, A, rhs, state.inverse_diagonal, T(0.7), rows)
            KP.gauss_seidel_sweep!(candidate, A, rhs, state.inverse_diagonal, T(0.7), rows)
            @test candidate == reference
            reference, candidate = copy(initial), copy(initial)
            invoke(KP.gauss_seidel_correction!, Tuple{Any, Any, Any, Any, Any, Any, Any},
                reference, state.correction, A, rhs, state.inverse_diagonal, T(0.7), rows)
            KP.gauss_seidel_correction!(candidate, state.correction, A, rhs,
                state.inverse_diagonal, T(0.7), rows)
            @test candidate == reference
        end
    end
end
