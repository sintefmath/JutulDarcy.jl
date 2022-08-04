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
            @test isapprox(p, p_ref, rtol = 1e-6)
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
            @test isapprox(p, p_ref, rtol = 1e-6)
            z = inj[:OverallMoleFractions]
            z_ref = repeat(i_mix, 1, 2)
            @test isapprox(z, z_ref, atol = 1e-8)
        end
    end
end
test_compositional_with_wells()
