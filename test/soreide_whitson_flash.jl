using JutulDarcy, MultiComponentFlash, StaticArrays, Test

@testset "SW flash reuse and phase identity" begin
    MCF = MultiComponentFlash
    mixture = MultiComponentMixture(["Water", "Methane", "Hydrogen"])
    eos = make_eos_immutable(GenericCubicEOS(mixture, SoreideWhitson(mixture)))
    full = JutulDarcy.FlashResults{SSIFlash, false, false}(SSIFlash(), false, 1e-8, 10.0, false, false)
    reuse = JutulDarcy.FlashResults{SSIFlash, false, true}(SSIFlash(), false, 1e-8, 10.0, false, true)
    initial = JutulDarcy.static_flashed_mixture(eos, Float64)
    cond = (p = 9e6, T = 320.0, z = @SVector [0.8, 0.18, 0.02])
    previous = JutulDarcy.immutable_flash_result(initial, full, eos, cond.p, cond.T, cond.z, 0.0)
    nearby = (p = 9.001e6, T = 320.01, z = @SVector [0.80001, 0.17999, 0.02])
    for options in (full, reuse)
        flash = JutulDarcy.immutable_flash_result(previous, options, eos,
            nearby.p, nearby.T, nearby.z, 0.0)
        @test flash.state == MCF.two_phase_lv
        @test flash.V ≈ 0.2 atol = 0.002
        x, y = flash.liquid.mole_fractions, flash.vapor.mole_fractions
        @test x[1] > 0.99
        @test y[1] < 0.01
        @test sum(x) ≈ 1 atol = 1e-10
        @test sum(y) ≈ 1 atol = 1e-10
        @test (1 - flash.V)*x + flash.V*y ≈ nearby.z atol = 1e-14
        forces = MCF.static_force_coefficients(eos, nearby, Float64)
        fl = MCF.static_fugacities(eos, (; nearby..., z = x, phase = Val(:liquid)), forces, Float64)
        fv = MCF.static_fugacities(eos, (; nearby..., z = y, phase = Val(:vapor)), forces, Float64)
        @test maximum(abs, log.(fl ./ fv)) < 1e-8
    end
    # Check sensitivities against fresh flashes, independently of the reuse
    # Newton Hessian and implicit derivative Jacobian.
    function vapor_fraction(u, options)
        p, T, water, methane = u
        z = SVector(water, methane, 1 - water - methane)
        return JutulDarcy.immutable_flash_result(previous, options, eos, p, T, z, 0.0).V
    end
    u = @SVector [nearby.p, nearby.T, nearby.z[1], nearby.z[2]]
    derivative = JutulDarcy.ForwardDiff.gradient(u -> vapor_fraction(u, reuse), u)
    for (i, h) in enumerate((100.0, 1e-3, 1e-5, 1e-5))
        delta = SVector{4}(ntuple(j -> i == j ? h : 0.0, 4))
        finite_difference = (vapor_fraction(u + delta, full) - vapor_fraction(u - delta, full))/(2h)
        @test derivative[i] ≈ finite_difference rtol = 1e-5 atol = 1e-10
    end
    for (z, state) in ((@SVector([0.001, 0.4995, 0.4995]), MCF.single_phase_v),
                      (@SVector([1 - 2e-10, 1e-10, 1e-10]), MCF.single_phase_l))
        flash = JutulDarcy.immutable_flash_result(previous, reuse, eos, 9e6, 333.15, z, 0.0)
        @test flash.state == state
        @test isfinite(flash.V)
    end
    # This stationary split was returned by the old common-reference flash.
    reversed_K = @SVector [673.7065250377914, 0.0001414624764567981, 6.720192989678175e-5]
    state, V, x, y, Z_l, Z_v = JutulDarcy.two_phase_flash_result(
        eos, cond, cond, 0.7998099354591234, reversed_K)
    reversed = MCF.FlashedMixture2Phase(state, reversed_K, V, x, y, Z_l, Z_v)
    recovered = JutulDarcy.immutable_flash_result(reversed, reuse, eos,
        cond.p + 100.0, cond.T, cond.z, 0.0)
    @test recovered.V ≈ previous.V atol = 1e-6
    @test recovered.liquid.mole_fractions[1] > 0.99
    @test recovered.vapor.mole_fractions[1] < 0.01
end
