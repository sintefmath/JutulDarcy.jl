using Jutul, JutulDarcy
using MultiComponentFlash
using Test

function test_compositional_with_wells()
    ## Define the mesh
    nx = 3
    ny = 1
    nz = 1
    # Some useful constants
    day = 3600*24
    bar = 1e5
    # Create and plot the mesh
    dims = (nx, ny, nz)
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    K = repeat([1e-13], number_of_cells(g))
    ## Set up a vertical well in the first corner, perforated in all layers
    prod = setup_vertical_well(g, K, nx, 1, name = :Producer);
    ## Set up an injector in the upper left corner
    inj = setup_vertical_well(g, K, 1, 1, name = :Injector);

    co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
    c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)
    c10 = MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488)

    mixture = MultiComponentMixture([co2, c1, c10], names = ["CO2", "C1", "C10"])
    eos = GenericCubicEOS(mixture)
    nc = number_of_cells(g)
    # Definition of fluid phases
    rhoLS, rhoVS = 1000.0, 100.0
    rhoS = [rhoLS, rhoVS]
    L, V = LiquidPhase(), VaporPhase()
    # Define system and realize on grid
    sys = MultiPhaseCompositionalSystemLV(eos, (L, V))
    model, parameters = setup_reservoir_model(g, sys, wells = [inj, prod], reference_densities = rhoS);
    state0 = setup_reservoir_state(model, Pressure = 150*bar, OverallMoleFractions = [0.5, 0.3, 0.2])

    dt = repeat([30.0]*day, 12*5)
    reservoir = reservoir_model(model);
    pv = pore_volume(model)
    inj_rate = 0.25*sum(pv)/sum(dt)
    rate_target = TotalRateTarget(inj_rate)
    i_mix =  [1.0, 0.0, 0.0]
    I_ctrl = InjectorControl(rate_target, i_mix, density = rhoVS)
    bhp_target = BottomHolePressureTarget(50*bar)
    P_ctrl = ProducerControl(bhp_target)

    controls = Dict()
    controls[:Injector] = I_ctrl
    controls[:Producer] = P_ctrl
    # Simulate
    forces = setup_reservoir_forces(model, control = controls)
    sim, config = setup_reservoir_simulator(model, state0, parameters, info_level = -1)
    states, reports = simulate!(sim, dt, forces = forces, config = config);

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
            z_ref = repeat(i_mix, 1, 2)
            @test isapprox(z, z_ref, atol = 1e-8)
        end
    end
end

function test_immiscible_with_wells()
    ## Define and plot the mesh
    nx = 3
    ny = 1
    nz = 1
    # Some useful constants
    day = 3600*24
    bar = 1e5
    # Create and plot the mesh
    dims = (nx, ny, nz)
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    ## Create a layered permeability field
    Darcy = 9.869232667160130e-13
    nlayer = nx*ny
    K = repeat([0.65*Darcy], nx*ny*nz)
    ## Set up a vertical well in the first corner, perforated in all layers
    P = setup_vertical_well(g, K, 1, 1, name = :Producer);
    ## Set up an injector in the upper left corner
    I = setup_well(g, K, [(nx, ny, 1)], name = :Injector);
    ## Set up a two-phase immiscible system and define a density secondary variable
    phases = (LiquidPhase(), VaporPhase())
    sys = ImmiscibleSystem(phases)
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    ## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
    model, parameters = setup_reservoir_model(g, sys, wells = [I, P], reference_densities = rhoS)
    ## Replace the density function with our custom version
    replace_variables!(model, PhaseMassDensities = ρ)
    ## Set up initial state
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
    ## Set up time-steps
    dt = repeat([30.0]*day, 12*5)
    pv = pore_volume(model)
    inj_rate = sum(pv)/sum(dt)
    rate_target = TotalRateTarget(inj_rate)
    i_mix = [0.0, 1.0]
    I_ctrl = InjectorControl(rate_target, i_mix, density = rhoGS)
    # The producer operates at a fixed bottom hole pressure
    bhp_target = BottomHolePressureTarget(50*bar)
    P_ctrl = ProducerControl(bhp_target)
    # Set up the controls. One control per well in the Facility.
    controls = Dict()
    controls[:Injector] = I_ctrl
    controls[:Producer] = P_ctrl
    # Set up forces for the whole model. For this example, all forces are defaulted
    # (amounting to no-flow for the reservoir).
    forces = setup_reservoir_forces(model, control = controls)
    ## Finally simulate!
    sim, config = setup_reservoir_simulator(model, state0, parameters, info_level = -1)
    states, reports = simulate!(sim, dt, forces = forces, config = config);

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
            s_ref = repeat(i_mix, 1, 2)
            @test isapprox(s, s_ref, atol = 1e-8)
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
    sys = ImmiscibleSystem(phases)
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    model, parameters = setup_reservoir_model(g, sys, wells = [I, P], reference_densities = rhoS)
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
test_compositional_with_wells()
test_immiscible_with_wells()
test_perforation_mask()
