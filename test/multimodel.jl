using Jutul, JutulDarcy
using Test

import JutulDarcy: simulate_mini_wellcase

function test_compositional_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:compositional_2ph_3c); kwarg...)
    @testset "Compositional with wells" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [5.341976350091003e6, 5.288647141719963e6, 5.221199515744034e6]
            @test isapprox(p, p_ref, rtol = 1e-4)
            z = res[:OverallMoleFractions]
            z_ref = [0.600928  0.52301   0.500993
                    0.145798  0.184558  0.193206
                    0.253274  0.292432  0.305802]
            @test isapprox(z, z_ref, atol = 1e-4)
        end

        @testset "Injector" begin
            inj = states[end][:Injector]
            p = inj[:Pressure]
            p_ref = [5.476226607917717e6, 5.476163345412835e6]
            @test isapprox(p, p_ref, rtol = 1e-4)
            z = inj[:OverallMoleFractions]
            z_ref = repeat([1.0, 0.0, 0.0], 1, 2)
            @test isapprox(z, z_ref, atol = 1e-8)
        end
    end
end

function test_immiscible_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:immiscible_2ph); kwarg...)
    @testset "Immiscible with wells" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [ 8.67141998033803e6, 1.7206847840192154e7, 2.567442267889495e7]
            @test isapprox(p, p_ref, rtol = 1e-4)
        end

        @testset "Injector" begin
            inj = states[end][:Injector]
            p = inj[:Pressure]
            p_ref = [ 2.9273662323437214e7, 2.9272083186770644e7]
            @test isapprox(p, p_ref, rtol = 1e-4)
            s = inj[:Saturations]
            s_ref = repeat([0.0, 1.0], 1, 2)
            @test isapprox(s, s_ref, atol = 1e-8)
        end
    end
end

function test_blackoil_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:bo_spe1); kwarg...)
    @testset "Blackoil with SPE1 PVT" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [5.05087e6, 5.14905e6, 5.21035e6]
            @test isapprox(p, p_ref, rtol = 1e-4)

            sw = res[:ImmiscibleSaturation]
            sw_ref = [0.0904386, 0.0884205, 0.0850772]
            @test isapprox(sw, sw_ref, atol = 1e-3)

            bo = states[end][:Reservoir][:BlackOilUnknown]
            for i in eachindex(bo)
                @test bo[i][2] == JutulDarcy.OilAndGas
            end
        end
    end
end

using JutulDarcy, Test
function test_perforation_mask()
    nx = 3
    ny = 1
    nz = 2
    day = 3600*24
    bar = 1e5
    dims = (nx, ny, nz)
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    Darcy = 9.869232667160130e-13
    K = repeat([0.65*Darcy], nx*ny*nz)
    P = setup_vertical_well(g, K, 1, 1, name = :Producer);
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    model, parameters = setup_reservoir_model(g, sys, wells = [P])
    replace_variables!(model, PhaseMassDensities = ρ)
    ## Set up initial state
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
    ## Set up time-steps
    dt = [30.0]*day
    # The producer operates at a fixed bottom hole pressure
    bhp_target = BottomHolePressureTarget(50*bar)
    P_ctrl = ProducerControl(bhp_target)
    # Set up the controls. One control per well in the Facility.
    controls = Dict()
    controls[:Producer] = P_ctrl
    forces = setup_reservoir_forces(model, control = controls)
    ## Mask away second perforation by multiplier
    pmask = PerforationMask([1.0, 0.0])
    forces[:Producer] = setup_forces(model.models[:Producer], mask = pmask)
    ## Simulate
    sim, config = setup_reservoir_simulator(model, state0, parameters, info_level = -1)
    states, reports = simulate!(sim, dt, forces = forces, config = config);
    v = states[1][:Producer][:TotalMassFlux]
    @testset "Perforation mask" begin
        @test abs(v[2]) < 1e-10
        @test abs(v[1]) > 1e-4
    end
end

@testset "MultiModel (wells)" begin 
    for gen_ad in [true, false]
        test_compositional_with_wells(general_ad = gen_ad)
        test_immiscible_with_wells(general_ad = gen_ad)
        test_blackoil_with_wells(general_ad = gen_ad)
    end
    test_perforation_mask()
end