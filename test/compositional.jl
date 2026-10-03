using JutulDarcy
using Jutul
using Test
using LinearAlgebra
using MultiComponentFlash
using StaticArrays

@testset "SIMPLE_COMP validation" begin
    dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("SIMPLE_COMP")
    data_path = joinpath(dpth, "SIMPLE_COMP.DATA")
    ref_path = joinpath(dpth, "reference.jld2")
    ref = Jutul.JLD2.load(ref_path)
    states_cmp = ref["e300"]
    for fast_flash in (false, true)
        case = setup_case_from_data_file(data_path, fast_flash = fast_flash)
        result = simulate_reservoir(case, info_level = -1)
        ws, states = result
        @testset "OverallMoleFractions" begin
            for i in 1:500
                rstate = states_cmp[i]
                state = states[i]
                z0 = rstate[:OverallMoleFractions]
                z = state[:OverallMoleFractions]
                @test norm(z-z0)/norm(z) < 0.05
            end
        end
        @testset "Pressure" begin
            for i in 1:500
                rstate = states_cmp[i]
                state = states[i]
                p0 = rstate[:Pressure]
                p = state[:Pressure]
                @test norm(p-p0)/norm(p) < 0.01
            end
        end

        @testset "Saturations" begin
            for i in 1:500
                rstate = states_cmp[i]
                state = states[i]
                s0 = rstate[:Saturations][2, :]
                s = state[:Saturations][2, :]
                @test norm(s-s0, 1)/norm(s, 1) < 0.05
            end
        end
    end
end

@testset "Independent compositional flash shortcuts" begin
    dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("SIMPLE_COMP")
    data_path = joinpath(dpth, "SIMPLE_COMP.DATA")

    # Check the model setup as well as the flash itself: neither option should
    # implicitly enable the other one.
    configurations = Dict{Tuple{Bool, Bool}, Any}()
    eos = nothing
    model = nothing
    for (reuse, bypass) in ((true, false), (false, true))
        case = setup_case_from_data_file(data_path;
            skip_wells = true,
            flash_reuse_guess = reuse,
            flash_stability_bypass = bypass)
        model = reservoir_model(case)
        flash = model.secondary_variables[:FlashResults]
        @test flash.reuse_guess == reuse
        @test flash.stability_bypass == bypass
        configurations[(reuse, bypass)] = flash
        eos = model.system.equation_of_state
    end
    reuse_only = configurations[(true, false)]
    bypass_only = configurations[(false, true)]
    both = JutulDarcy.FlashResults(model;
        reuse_guess = true, stability_bypass = true)
    z = @SVector [0.6, 0.1, 0.3]
    initial = JutulDarcy.static_flashed_mixture(eos, Float64)

    previous_two_phase = JutulDarcy.immutable_flash_result(
        initial, reuse_only, eos, 1.0e6, 300.0, z, 0.0)
    @test previous_two_phase.state == MultiComponentFlash.two_phase_lv
    V0 = previous_two_phase.V
    @test JutulDarcy.estimated_K_from_previous_flash(
        previous_two_phase, V0, z) ≈ previous_two_phase.K rtol = 1e-8
    z_near = @SVector [0.599, 0.101, 0.3]
    estimated_K = JutulDarcy.estimated_K_from_previous_flash(
        previous_two_phase, V0, z_near)
    @test all(isfinite, estimated_K)
    @test norm(estimated_K - previous_two_phase.K) > 1e-5

    nearby_two_phase = (p = 1.01e6, T = 300.0, z = z_near)
    V_reused, K_reused, reused_stability = JutulDarcy.numeric_flash(
        previous_two_phase, reuse_only, eos, nearby_two_phase)
    V_full, K_full, full_stability = JutulDarcy.numeric_flash(
        previous_two_phase, bypass_only, eos, nearby_two_phase)
    @test 0.0 < V_reused < 1.0
    @test V_reused ≈ V_full rtol = 1e-7
    @test K_reused ≈ K_full rtol = 1e-7
    @test !reused_stability.report.liquid.trivial # quick two-phase attempt
    @test full_stability.report.liquid.trivial # full stability test
    @test !full_stability.bypassed # no valid single-phase bypass cache

    # A phase change must reject the quick result and use a full flash.
    single_phase = (p = 7.5e6, T = 300.0, z = z)
    V_transition, _, transition_stability = JutulDarcy.numeric_flash(
        previous_two_phase, reuse_only, eos, single_phase)
    @test isnan(V_transition)
    @test transition_stability.stable

    previous_single_phase = JutulDarcy.immutable_flash_result(
        initial, bypass_only, eos, 1.0e5, 800.0, z, 0.0)
    @test previous_single_phase.state == MultiComponentFlash.single_phase_v
    @test isfinite(previous_single_phase.critical_distance)
    near_stable = (p = 1.001e5, T = 800.01,
        z = @SVector [0.60001, 0.09999, 0.3])
    _, _, bypassed = JutulDarcy.numeric_flash(
        previous_single_phase, bypass_only, eos, near_stable)
    _, _, tested = JutulDarcy.numeric_flash(
        previous_single_phase, reuse_only, eos, near_stable)
    @test bypassed.bypassed
    @test bypassed.storage.reference == previous_single_phase.flash_cond
    @test !tested.bypassed
    @test isnan(tested.storage.critical_distance)

    # With both shortcuts on, the prior phase state still determines which one
    # applies; they do not depend on one another.
    _, _, combined_two_phase = JutulDarcy.numeric_flash(
        previous_two_phase, both, eos, nearby_two_phase)
    _, _, combined_single_phase = JutulDarcy.numeric_flash(
        previous_single_phase, both, eos, near_stable)
    @test !combined_two_phase.bypassed
    @test !combined_two_phase.report.liquid.trivial
    @test combined_single_phase.bypassed
end
