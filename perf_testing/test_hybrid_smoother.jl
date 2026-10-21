using Jutul, Test, SparseArrays, LinearAlgebra, Random
const KP = Jutul.KAPreconditioners
Base.include(KP, joinpath(@__DIR__, "hybrid_smoother_candidate.jl"))

@testset "Hybrid uses existing level smoothers" begin
    rng = MersenneTwister(41)
    for T in (Float32, Float64), steps in (1, 2)
        n = 4096
        M = spdiagm(-64 => fill(T(-0.8), n - 64), -1 => fill(T(-0.9), n - 1),
            0 => fill(T(5), n), 1 => fill(T(-1), n - 1), 64 => fill(T(-0.7), n - 64))
        A = KP.csr_matrix(M)
        config = KP.PerfHybridSmoother(1000; steps)
        options = KP.AMGOptions(; smoother = config)
        H = KP.setup_amg(A, options)
        reference = KP.setup_amg(A, KP.AMGOptions(; smoother = KP.SPAI0(steps)))
        @test length(H.levels) == length(reference.levels)
        for (level, original) in zip(H.levels, reference.levels)
            fine = size(level.A, 1) > config.cutoff
            @test fine ? level.smoother isa KP.SPAI0State : level.smoother isa KP.GaussSeidelState
            @test size(level.A) == size(original.A)
            if !fine
                original.smoother = KP.setup_smoother(original.A, KP.GaussSeidel(steps))
            end
        end
        rhs = randn(rng, T, n)
        x, expected = similar(rhs), similar(rhs)
        KP.apply!(expected, reference, rhs)
        KP.apply!(x, H, rhs)
        @test x == expected
        diagonal = H.levels[1].smoother.diagonal
        B = KP.csr_matrix(copy(M))
        B.nzval .*= T(2)
        KP.resetup_amg!(H, B, :operators)
        KP.apply!(x, H, rhs)
        @test x ≈ expected/T(2)
        @test H.levels[1].smoother.diagonal === diagonal
        KP.resetup_amg!(H, B, :memory)
        KP.apply!(x, H, rhs)
        @test all(isfinite, x)
        @test H.levels[1].smoother.diagonal === diagonal
    end
end

include("test_smoother_reuse.jl")
