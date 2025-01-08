using Jutul, JutulDarcy
using Test

import JutulDarcy: simulate_mini_wellcase

function test_compositional_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:compositional_2ph_3c); simple_well = false, kwarg...)
    @testset "Compositional with wells" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [5.218139052839837e6, 5.180811416153696e6, 5.137566116526595e6]
            @test isapprox(p, p_ref, rtol = 1e-3)
            z = res[:OverallMoleFractions]
            z_ref = [
                0.6032429704979841 0.521568527073734 0.5000266982944569;
                0.13981926168012976 0.181469203945046 0.1907407040733318; 
                0.2569377678218861 0.2969622689812199 0.3092325976322113
                ]
            @test isapprox(z, z_ref, atol = 1e-3)
        end

        @testset "Injector" begin
            inj = states[end][:Injector]
            p = inj[:Pressure]
            p_ref = [5.332306160997921e6, 5.332306160997921e6]
            @test isapprox(p, p_ref, rtol = 1e-3)
            z = inj[:OverallMoleFractions]
            z_ref = repeat([1.0, 0.0, 0.0], 1, 2)
            @test isapprox(z, z_ref, atol = 1e-8)
        end
    end
end

function test_immiscible_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:immiscible_2ph); simple_well = false, kwarg...)
    @testset "Immiscible with wells" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [2.8349045684888966e7, 3.67188418184439e7, 4.50238389828125e7]
            @test isapprox(p, p_ref, rtol = 1e-4)
        end

        @testset "Injector" begin
            inj = states[end][:Injector]
            p = inj[:Pressure]
            p_ref = [6.753641685512246e7, 6.753793670867026e7]
            @test isapprox(p, p_ref, rtol = 1e-4)
            s = inj[:Saturations]
            s_ref = repeat([0.0, 1.0], 1, 2)
            @test isapprox(s, s_ref, atol = 1e-8)
        end
    end
end

function test_blackoil_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:bo_spe1); simple_well = false, kwarg...)
    @testset "Blackoil with SPE1 PVT" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [5.305212867881992e6, 5.409995957652735e6, 5.472027929725644e6]
            @test isapprox(p, p_ref, rtol = 1e-4)

            sw = res[:ImmiscibleSaturation]
            sw_ref = [0.09318952675271851, 0.09143892412276486, 0.08736801209544151]
            @test isapprox(sw, sw_ref, atol = 1e-3)

            bo = states[end][:Reservoir][:BlackOilUnknown]
            for i in eachindex(bo)
                @test bo[i].phases_present == JutulDarcy.OilAndGas
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
    domain = reservoir_domain(g, permeability = 0.1*Darcy, porosity = 0.1)
    P = setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = false);
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    visLS = 1e-4
    visGS = 1e-3
    parameters = Dict(:Reservoir=>Dict(:PhaseViscosities=>[visLS, visGS]))
    model, parameters = setup_reservoir_model(domain, sys, wells = [P], parameters=parameters)
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
        @test abs(v[2]) < 1e-12
        @test abs(v[1]) > 1e-7
    end
end

@testset "MultiModel (wells)" begin
    for b in [false, true]
        for backend in [:csr, :csc]
            for gen_ad in [true, false]
                @testset "Block=$b, backend=$b" begin
                    arg = (general_ad = gen_ad, backend = backend, block_backend = b)
                    test_compositional_with_wells(; arg...)
                    test_compositional_with_wells(; fast_flash = true, arg...)
                    test_immiscible_with_wells(; arg...)
                    test_blackoil_with_wells(; arg...)
                end
            end
        end
    end
    test_perforation_mask()

    @testset "Preconditioners" begin
        # CSR is the default, do exhaustive testing
        arg = (
            general_ad = false,
            backend = :csr,
            block_backend = true
        )
        for lsolve in [:bicgstab, :gmres]
            for prec in [:cpr, :cprw, :ilu0, :jacobi, :spai0]
                @testset "$lsolve:$prec" begin
                    lsolve_arg = (
                        precond = prec,
                        linear_solver = lsolve,
                    )
                    test_compositional_with_wells(; setuparg = lsolve_arg, arg...)
                    test_immiscible_with_wells(; setuparg = lsolve_arg, arg...)
                    test_blackoil_with_wells(; setuparg = lsolve_arg, arg...)
                end
            end
        end
    end
end
