using Jutul, JutulDarcy
using Test

import JutulDarcy: simulate_mini_wellcase

function test_compositional_with_wells(; kwarg...)
    states, = simulate_mini_wellcase(Val(:compositional_2ph_3c); simple_well = false, kwarg...)
    @testset "Compositional with wells" begin
        @testset "Reservoir" begin
            res = states[end][:Reservoir]
            p = res[:Pressure]
            p_ref = [5.218139052839806e6, 5.180811416153673e6, 5.137566116526581e6]
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
            p_ref = [5.333171412614246e6, 5.333171412614246e6]
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
            p_ref = [2.834901722247195e7, 3.6718809754436985e7, 4.5023807041482836e7]
            @test isapprox(p, p_ref, rtol = 1e-4)
        end

        @testset "Injector" begin
            inj = states[end][:Injector]
            p = inj[:Pressure]
            p_ref = [6.753082589164461e7, 6.753082589164461e7]
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
            p_ref = [5.305212863516781e6, 5.409996331019564e6, 5.472028582373638e6]
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
    P = setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = false, use_top_node = true);
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
    model, parameters = setup_reservoir_model(domain, sys, wells = [P], parameters=parameters, extra_out = true)
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

function test_group_control()
    day = 3600*24
    bar = 1e5
    Darcy = 9.869232667160130e-13
    nx = ny = 5
    nz = 2
    g = CartesianMesh((nx, ny, nz), (2000.0, 1500.0, 50.0))
    domain = reservoir_domain(g, permeability = 0.3*Darcy, porosity = 0.2)
    Prod1 = setup_vertical_well(domain, 1, 1, name = :Producer1)
    Prod2 = setup_vertical_well(domain, nx, ny, name = :Producer2)
    Inj = setup_well(domain, [(nx÷2+1, ny÷2+1, 1)], name = :Injector)
    groups = Dict(:ProducerGroup => [:Producer1, :Producer2])
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    model, parameters = setup_reservoir_model(domain, sys,
        wells = [Inj, Prod1, Prod2],
        extra_out = true,
        well_groups = groups
    )
    replace_variables!(model, PhaseMassDensities = ρ)
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
    dt = repeat([30.0]*day, 6)
    pv = pore_volume(model, parameters)
    inj_rate = 0.25*sum(pv)/sum(dt)
    rate_target = TotalRateTarget(inj_rate)
    I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
    # Producers use GroupTarget - the actual target comes from group_controls
    P1_ctrl = ProducerControl(GroupTarget(:ProducerGroup))
    P2_ctrl = ProducerControl(GroupTarget(:ProducerGroup))
    controls = Dict{Symbol, Any}()
    controls[:Injector] = I_ctrl
    controls[:Producer1] = P1_ctrl
    controls[:Producer2] = P2_ctrl
    # Group controls define the actual target for the group
    group_controls = Dict{Symbol, Any}(
        :ProducerGroup => ProducerControl(TotalRateTarget(-inj_rate))
    )
    forces = setup_reservoir_forces(model, control = controls, group_controls = group_controls)
    result = simulate_reservoir(state0, model, dt,
        parameters = parameters, forces = forces, info_level = -1
    )
    wd, states, t = result
    @testset "Group control simulation" begin
        @test length(states) == 6
        # Both producers should produce approximately equal rates under group control
        q1 = wd[:Producer1][:rate]
        q2 = wd[:Producer2][:rate]
        for step in eachindex(q1)
            if q1[step] != 0.0 && q2[step] != 0.0
                ratio = q1[step] / q2[step]
                @test isapprox(ratio, 1.0, atol = 0.15)
            end
        end
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

    @testset "Group control" begin
        test_group_control()
    end

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
