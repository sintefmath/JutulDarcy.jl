function simulate_mini_wellcase(::Val{:compositional_2ph_3c};
        dims = (3, 1, 1),
        setuparg = NamedTuple(),
        output_path = nothing,
        default_linsolve = true,
        nstep = 12*5,
        total_time = 30.0*si_unit(:day)*nstep,
        kwarg...
    )
    # Some useful constants
    day = 3600*24
    bar = 1e5
    # Create the mesh
    nx, ny, nz = dims
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    domain = reservoir_domain(g, permeability = 1e-13, porosity = 0.1)
    ## Set up a vertical well in the first corner, perforated in all layers
    prod = setup_vertical_well(domain, nx, 1, name = :Producer);
    ## Set up an injector in the upper left corner
    inj = setup_vertical_well(domain, 1, 1, name = :Injector);

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
    model, parameters = setup_reservoir_model(domain, sys, wells = [inj, prod]; kwarg...);
    state0 = setup_reservoir_state(model, Pressure = 150*bar, OverallMoleFractions = [0.5, 0.3, 0.2])

    dt = fill(total_time/nstep, nstep)
    pv = pore_volume(domain)
    time_scale = 30.0*12*5*si_unit(:day)
    inj_rate = 0.25*sum(pv)/time_scale
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
    sim, config = setup_reservoir_simulator(
        model, state0, parameters;
        info_level = -1,
        set_linear_solver = !default_linsolve,
        error_on_incomplete = true,
        output_path = output_path,
        setuparg...,
        tol_cnv_well = 1e-3 # Needed for test
    )
    setup = Dict(:config     => config,
                 :forces     => forces,
                 :state0     => state0,
                 :model      => model,
                 :sim        => sim,
                 :parameters => parameters,
                 :dt         => dt)
    states, reports = simulate!(sim, dt, forces = forces, config = config);
    return (states = states, reports = reports, setup = setup)
end

function simulate_mini_wellcase(::Val{:immiscible_2ph};
        dims = (3, 1, 1),
        setuparg = NamedTuple(),
        output_path = nothing,
        permeability = 0.1*9.869232667160130e-13,
        nstep = 12*5,
        total_time = 30.0*si_unit(:day)*nstep,
        default_linsolve = true,
        kwarg...)
    # Some useful constants
    day = 3600*24
    bar = 1e5
    # Create the mesh
    nx, ny, nz = dims
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    domain = reservoir_domain(g, permeability = permeability, porosity = 0.1)
    ## Set up a vertical well in the first corner, perforated in all layers
    P = setup_vertical_well(domain, 1, 1, name = :Producer);
    ## Set up an injector in the upper left corner
    I = setup_well(domain, [(nx, ny, 1)], name = :Injector);
    ## Set up a two-phase immiscible system and define a density secondary variable
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    ## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
    model, parameters = setup_reservoir_model(domain, sys, wells = [I, P]; kwarg...)
    ## Replace the density function with our custom version
    replace_variables!(model, PhaseMassDensities = ρ)
    ## Set up initial state
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
    ## Set up time-steps
    dt = fill(total_time/nstep, nstep)
    pv = pore_volume(model, parameters)
    time_scale = 30.0*12*5*si_unit(:day)
    inj_rate = sum(pv)/time_scale
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
    sim, config = setup_reservoir_simulator(
        model, state0, parameters;
        info_level = -1,
        set_linear_solver = !default_linsolve,
        output_path = output_path,
        error_on_incomplete = true,
        setuparg...
        )
    setup = Dict(:config     => config,
                 :forces     => forces,
                 :state0     => state0,
                 :model      => model,
                 :sim        => sim,
                 :parameters => parameters,
                 :dt         => dt)
    states, reports = simulate!(sim, dt, forces = forces, config = config);
    return (states = states, reports = reports, setup = setup)
end

function simulate_mini_wellcase(::Val{:bo_spe1};
        dims = (3, 1, 1),
        setuparg = NamedTuple(),
        output_path = nothing,
        default_linsolve = true,
        nstep = 12*5,
        total_time = 30.0*si_unit(:day)*nstep,
        kwarg...
    )
    # Some useful constants
    day = 3600*24
    bar = 1e5
    # Create the mesh
    nx, ny, nz = dims
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    ## Create a layered permeability field
    Darcy = 9.869232667160130e-13
    domain = reservoir_domain(g, permeability = 0.1*Darcy, porosity = 0.1)
    ## Set up a vertical well in the first corner, perforated in all layers
    P = setup_vertical_well(domain, 1, 1, name = :Producer);
    ## Set up an injector in the upper left corner
    I = setup_well(domain, [(nx, ny, 1)], name = :Injector);
    ## Set up a two-phase immiscible system and define a density secondary variable
    setup = JutulDarcy.blackoil_bench_pvt(:spe1)
    pvt = setup[:pvt]
    pvto = pvt[2]
    rhoS = setup[:rhoS]
    phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
    pvto_tab = pvto.tab[1]
    sat_table = get_1d_interpolator(pvto_tab.sat_pressure, pvto_tab.rs, cap_end = false)
    sys = StandardBlackOilSystem(rs_max = sat_table, phases = phases, reference_densities = rhoS)
    ## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
    model, parameters = setup_reservoir_model(domain, sys, wells = [I, P]; kwarg...)
    ## Set up initial state
    bo = BlackOilX(50.0, JutulDarcy.OilOnly, false)
    state0 = setup_reservoir_state(model, Pressure = 200*bar, ImmiscibleSaturation = 0.1, BlackOilUnknown = bo)
    ## Set up time-steps
    dt = fill(total_time/nstep, nstep)
    pv = pore_volume(domain)
    time_scale = 30.0*12*5*si_unit(:day)
    inj_rate = sum(pv)/time_scale
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
    sim, config = setup_reservoir_simulator(
        model, state0, parameters;
        info_level = -1,
        set_linear_solver = !default_linsolve,
        output_path = output_path,
        error_on_incomplete = true,
        setuparg...
        )
    setup = Dict(:config     => config,
                 :forces     => forces,
                 :state0     => state0,
                 :model      => model,
                 :sim        => sim,
                 :parameters => parameters,
                 :dt         => dt)
    states, reports = simulate!(sim, dt, forces = forces, config = config);
    return (states = states, reports = reports, setup = setup)
end

function simulate_mini_wellcase(physics::Symbol; kwarg...)
    return simulate_mini_wellcase(Val(physics); kwarg...)
end

function precompile_darcy_multimodels(targets = missing; kwarg...)
    if ismissing(targets)
        targets = []
        # Block backend, CSC (CPR)
        push!(targets, (true, :csc))
        # Scalar blackend, CSC (direct solver)
        # push!(targets, (false, :csc))
        # Block backend, CSR (CPR)
        # push!(targets, (true, :csr))
    end
    wellcases = [:bo_spe1, :immiscible_2ph, :compositional_2ph_3c]
    for (block_backend, backend) in targets
        for w in wellcases
            simulate_mini_wellcase(Val(w); block_backend = block_backend, backend = backend, kwarg...)
        end
    end
    try
        simulate_mini_wellcase(Val(:immiscible_2ph); output_path = tempdir()) 
    catch
        # Don't do anything, cannot write to folder?
    end
    return true
end
