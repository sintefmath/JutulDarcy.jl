using Jutul, JutulDarcy, MultiComponentFlash, StaticArrays, Test
using KernelAbstractions, JLArrays

@kernel function phase_amount_reuse_kernel!(out, derivative, old, settings, eos, cond)
    i = @index(Global)
    p = JutulDarcy.ForwardDiff.Dual{Nothing}(cond.p, one(cond.p))
    result = JutulDarcy.immutable_flash_result(old, settings, eos, p, cond.T, cond.z, 0.0)
    @inbounds out[i] = JutulDarcy.compositional_primal(result.V)
    @inbounds derivative[i] = JutulDarcy.ForwardDiff.partials(result.V)[1]
end

@testset "Phase-amount Newton reuse" begin
    components = [
        MolecularProperty(0.044, 7.38e6, 304.1, 9.412e-5, 0.224),
        MolecularProperty(0.016, 4.6e6, 190.6, 9.863e-5, 0.011),
        MolecularProperty(0.0142, 2.1e6, 617.7, 6.098e-4, 0.488),
    ]
    eos = make_eos_immutable(GenericCubicEOS(MultiComponentMixture(components)))
    fresh = JutulDarcy.FlashResults{SSIFlash, false, false}(
        SSIFlash(), false, 1.0e-8, 10.0, false, false
    )
    reuse = JutulDarcy.FlashResults{SSIFlash, false, true}(
        SSIFlash(), false, 1.0e-8, 10.0, false, true
    )
    initial = JutulDarcy.static_flashed_mixture(eos, Float64)
    z0 = @SVector [0.5, 0.3, 0.2]
    old = JutulDarcy.immutable_flash_result(initial, fresh, eos, 1.0e6, 300.0, z0, 0.0)
    @test old.state == MultiComponentFlash.two_phase_lv
    z = @SVector [0.50001, 0.29999, 0.2]
    cond = (p = 1.0001e6, T = 300.01, z = z)
    liquid, vapor = JutulDarcy.estimated_phase_amounts_from_previous_flash(old, old.V, z)
    @test liquid + vapor ≈ z
    @test vapor ./ z ≈ old.V * old.vapor.mole_fractions ./ z0
    @test sum(liquid / sum(liquid)) ≈ 1
    @test sum(vapor / sum(vapor)) ≈ 1

    K0 = (vapor / sum(vapor)) ./ (liquid / sum(liquid))
    V0 = sum(vapor)
    V, K, report = flash_2ph!(
        MultiComponentFlash.StaticConfig(), K0, eos, cond, V0;
        method = NewtonFlash(dMax = 1.0), maxiter = 3, extra_out = true,
        tolerance = 1.0e-8, z_min = nothing, check = false
    )
    @test report.converged
    @test report.its <= 3
    actual = JutulDarcy.immutable_flash_result(old, reuse, eos, cond.p, cond.T, z, 0.0)
    expected = JutulDarcy.immutable_flash_result(old, fresh, eos, cond.p, cond.T, z, 0.0)
    @test actual.state == expected.state == MultiComponentFlash.two_phase_lv
    @test actual.V ≈ V atol = 1.0e-10
    @test actual.K ≈ K rtol = 1.0e-9
    @test actual.V ≈ expected.V atol = 1.0e-8
    @test actual.liquid.mole_fractions ≈ expected.liquid.mole_fractions rtol = 1.0e-7
    @test actual.vapor.mole_fractions ≈ expected.vapor.mole_fractions rtol = 1.0e-7

    if isdefined(JLArrays, :JLBackend)
        out = JLArray(zeros(8))
        derivative = JLArray(zeros(8))
        backend = JLArrays.JLBackend()
        phase_amount_reuse_kernel!(backend, 8)(
            out, derivative, old, reuse, eos, cond;
            ndrange = 8
        )
        p_ad = JutulDarcy.ForwardDiff.Dual{Nothing}(cond.p, 1.0)
        reference_ad = JutulDarcy.immutable_flash_result(
            old, reuse, eos, p_ad, cond.T, z, 0.0
        )
        @test all(isapprox.(Array(out), actual.V; atol = 1.0e-10, rtol = 0))
        @test all(isapprox.(Array(derivative), reference_ad.V.partials[1]; rtol = 1.0e-10))
    end

    # Leaving the two-phase region must use the full stability/SSI fallback.
    single = JutulDarcy.immutable_flash_result(old, reuse, eos, 1.0e5, 800.0, z, 0.0)
    single_reference = JutulDarcy.immutable_flash_result(old, fresh, eos, 1.0e5, 800.0, z, 0.0)
    @test single.state == single_reference.state
    @test single.V == single_reference.V

    # An extreme trace-component K does not make coincident phases reusable.
    @test !JutulDarcy.valid_reused_flash(
        0.3, @SVector([1.0, 1.0e28]), true,
        @SVector([1 - 1.0e-10, 1.0e-10]), 1.0e-8
    )
    @test !JutulDarcy.valid_reused_flash(V, K, false, z, 1.0e-8)

    # Do not replace a physically tiny old liquid amount by MINIMUM_COMPOSITION:
    # that changes theta_i when the feed later gains that component.
    trace_z = @SVector [0.8, 0.2 - 1.0e-10, 1.0e-10]
    trace_K = @SVector [0.1, 10.0, 1.0e28]
    trace_V = solve_rachford_rice(trace_K, trace_z)
    trace_x, trace_y = JutulDarcy.phase_mole_fractions(trace_z, trace_K, trace_V)
    trace = FlashedMixture2Phase(MultiComponentFlash.two_phase_lv, trace_K, trace_V,
        trace_x, trace_y, 1.0, 1.0, NaN, (; cond..., z = trace_z))
    trace_feed = @SVector [0.1, 0.89, 0.01]
    trace_liquid, trace_vapor = JutulDarcy.estimated_phase_amounts_from_previous_flash(
        trace, trace_V, trace_feed)
    @test trace_liquid + trace_vapor ≈ trace_feed
    @test 0 < trace_liquid[3] < 1.0e-26
    @test trace_vapor[3] ≈ trace_feed[3]
end
