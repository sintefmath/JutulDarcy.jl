using Test, JutulDarcy

@testset "Reservoir voidage rate derivatives" begin
    JD = JutulDarcy
    FD = JD.ForwardDiff
    bW, bO, bG = 1.02, 0.8, 0.01
    state_avg = (bW = bW, bO = rs -> bO, bG = rv -> bG, p = 100e5, rs = 100.0, rv = 1e-4)
    resv(q) = JD.compute_total_resv_rate(state_avg; qw = q[1], qo = q[2], qg = q[3])
    # A reopened producer starts at (almost) zero rates. The derivatives must
    # not depend on round-off in the rates.
    for q in ([0.0, 0.0, 0.0], [-1e-30, 0.0, 0.0], [0.0, -1e-30, 0.0])
        @test FD.gradient(resv, q) ≈ [1/bW, 1/bO, 1/bG]
    end
    # Producer and injector give the same volume with opposite sign
    q = [-1.0, -2.0, -300.0]
    @test resv(q) ≈ -resv(-q)
    @test resv(q) < 0
end
