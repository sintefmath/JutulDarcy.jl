function simulate_mini_wellcase(::Val{:compositional_2ph_3c}; kwarg...)
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
    # Definition of fluid phases
    rhoLS, rhoVS = 1000.0, 100.0
    rhoS = [rhoLS, rhoVS]
    L, V = LiquidPhase(), VaporPhase()
    # Define system and realize on grid
    sys = MultiPhaseCompositionalSystemLV(eos, (L, V), reference_densities = rhoS)
    model, parameters = setup_reservoir_model(g, sys, wells = [inj, prod]; kwarg...);
    state0 = setup_reservoir_state(model, Pressure = 150*bar, OverallMoleFractions = [0.5, 0.3, 0.2])

    dt = repeat([30.0]*day, 12*5)
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
    return simulate!(sim, dt, forces = forces, config = config);
end

function simulate_mini_wellcase(::Val{:immiscible_2ph}; kwarg...)
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
    K = repeat([0.65*Darcy], nx*ny*nz)
    ## Set up a vertical well in the first corner, perforated in all layers
    P = setup_vertical_well(g, K, 1, 1, name = :Producer);
    ## Set up an injector in the upper left corner
    I = setup_well(g, K, [(nx, ny, 1)], name = :Injector);
    ## Set up a two-phase immiscible system and define a density secondary variable
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    ## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
    model, parameters = setup_reservoir_model(g, sys, wells = [I, P]; kwarg...)
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
    return simulate!(sim, dt, forces = forces, config = config);
end

function simulate_mini_wellcase(::Val{:bo_spe1}; kwarg...)
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
    setup = JutulDarcy.blackoil_bench_pvt(:spe1)
    pvt = setup[:pvt]
    pvto = pvt[2]
    rhoS = setup[:rhoS]
    phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
    sat_table = get_1d_interpolator(pvto.sat_pressure, pvto.rs, cap_end = false)
    sys = StandardBlackOilSystem(sat_table, phases = phases, reference_densities = rhoS)
    ## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
    model, parameters = setup_reservoir_model(g, sys, wells = [I, P]; kwarg...)
    ## Set up initial state
    bo = (50.0, JutulDarcy.OilOnly, false)
    state0 = setup_reservoir_state(model, Pressure = 200*bar, ImmiscibleSaturation = 0.1, BlackOilUnknown = bo)
    ## Set up time-steps
    dt = repeat([30.0]*day, 12*5)
    pv = pore_volume(model)
    inj_rate = sum(pv)/sum(dt)
    rate_target = TotalRateTarget(inj_rate)
    i_mix = [0.01, 0.05, 0.94]
    I_ctrl = InjectorControl(rate_target, i_mix, density = 1.0)
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
    return simulate!(sim, dt, forces = forces, config = config);
end

function precompile_darcy_multimodels()
    targets = []
    # Block backend, CSC (CPR)
    push!(targets, (true, :csc))
    # Scalar blackend, CSC (direct solver)
    push!(targets, (false, :csc))
    # Block backend, CSR (CPR)
    # push!(targets, (true, :csr))
    wellcases = [:bo_spe1, :immiscible_2ph, :compositional_2ph_3c]
    for (block_backend, backend) in targets
        for w in wellcases
            simulate_mini_wellcase(Val(w), block_backend = block_backend, backend = backend)
        end
    end
end
