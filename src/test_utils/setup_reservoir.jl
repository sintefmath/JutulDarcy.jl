function get_test_setup(mesh_or_casename; case_name = "single_phase_simple", context = "cpu", timesteps = [1.0, 2.0], pvfrac = 0.05, kwarg...)
    if isa(mesh_or_casename, String)
        G, mrst_data = get_minimal_tpfa_grid_from_mrst(mesh_or_casename, extraout = true)
        mesh = MRSTWrapMesh(mrst_data["G"])
    else
        mesh = mesh_or_casename
        geo = tpfv_geometry(mesh)
        G = discretized_domain_tpfv_flow(geo; kwarg...)
    end
    nc = number_of_cells(G)
    pv = G.grid.pore_volumes
    timesteps = timesteps*3600*24

    if context == "cpu"
        context = DefaultContext()
    elseif isa(context, String)
        error("Unsupported target $context")
    end
    @assert isa(context, JutulContext)
    parameters = nothing
    if case_name == "single_phase_simple"
        # Parameters
        bar = 1e5
        p0 = 100*bar # 100 bar
        cl = 1e-5/bar
        pRef = 100*bar
        rhoLS = 1000
        # Single-phase liquid system (compressible pressure equation)
        phase = LiquidPhase()
        sys = SinglePhaseSystem(phase)
        # Simulation model wraps grid and system together with context (which will be used for GPU etc)
        model = SimulationModel(G, sys, context = context)
        s = model.secondary_variables
        s[:PhaseMassDensities] = ConstantCompressibilityDensities(sys, pRef, rhoLS, cl)

        # System state
        tot_time = sum(timesteps)
        irate = pvfrac*sum(pv)/tot_time
        src = SourceTerm(1, irate)
        bc = FlowBoundaryCondition(nc, p0/2)
        forces = setup_forces(model, sources = src, bc = bc)

        # State is dict with pressure in each cell
        init = Dict(:Pressure => p0)
    elseif case_name == "single_phase_fakewells"
        # Parameters
        bar = 1e5
        cl = 1e-5/bar
        pRef = 100*bar
        rhoLS = 1000
        # Single-phase liquid system (compressible pressure equation)
        phase = LiquidPhase()
        sys = SinglePhaseSystem(phase)
        # Simulation model wraps grid and system together with context (which will be used for GPU etc)
        model = SimulationModel(G, sys, context = context)
        s = model.secondary_variables
        s[:PhaseMassDensities] = ConstantCompressibilityDensities(sys, pRef, rhoLS, cl)
        forces = setup_forces(model)

        pv = model.domain.grid.pore_volumes
        pv[1] *= 1000
        pv[end] *= 1000

        p_init = 100*bar
        nc = number_of_cells(G)
        p0 = repeat([p_init], nc)
        p0[1] = 2*p_init
        p0[end] = 0.5*p_init
        init = Dict(:Pressure => p0)
    elseif case_name == "two_phase_simple"
        bar = 1e5
        p0 = 100*bar # 100 bar
        mu = 1e-3    # 1 cP
        cl = 1e-5/bar
        pRef = 100*bar
        rhoLS = 1000
        L = LiquidPhase()
        V = VaporPhase()
        sys = ImmiscibleSystem([L, V])
        model = SimulationModel(G, sys, context = context)

        kr = BrooksCoreyRelPerm(sys, [2, 3])
        s = model.secondary_variables
        s[:RelativePermeabilities] = kr
        s[:PhaseMassDensities] = ConstantCompressibilityDensities(sys, pRef, rhoLS, cl)
        parameters = setup_parameters(model, PhaseViscosities = [mu, mu/2])

        tot_time = sum(timesteps)
        irate = pvfrac*sum(pv)/tot_time
        s0 = 1.0
        s = 0.0

        # s = 0.1
        # s0 = 0.9
        src  = [SourceTerm(1, irate, fractional_flow = [1 - s, s]), 
                SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
        forces = setup_forces(model, sources = src)

        # State is dict with pressure in each cell
        init = Dict(:Pressure => p0, :Saturations => [1 - s0, s0])
    elseif case_name == "two_phase_fake_wells"
        inj = 1
        prod = nc
        G.grid.pore_volumes[inj] *= 1000
        G.grid.pore_volumes[prod] *= 1000

        bar = 1e5
        p0 = 100*bar # 100 bar
        s0 = 1.0

        mu = 1e-3    # 1 cP
        cl = 1e-5/bar
        pRef = 100*bar
        rhoLS = 1000
        L = LiquidPhase()
        V = VaporPhase()
        sys = ImmiscibleSystem([L, V])
        model = SimulationModel(G, sys, context = context)

        kr = BrooksCoreyRelPerm(sys, [2, 3])
        s = model.secondary_variables
        s[:RelativePermeabilities] = kr
        s[:PhaseMassDensities] = ConstantCompressibilityDensities(sys, pRef, rhoLS, cl)
        parameters = setup_parameters(model, PhaseViscosities = [mu, mu/2])

        tot_time = sum(timesteps)
        forces = setup_forces(model)

        p_init = repeat([p0], nc)
        p_init[inj] = 2*p0
        p_init[prod] = p0/2

        s_init = repeat([1 - s0, s0], 1, nc)
        s_init[1, inj] = s0
        s_init[2, inj] = 1 - s0

        # State is dict with pressure in each cell
        init = Dict(:Pressure => p_init, :Saturations => s_init)
    elseif case_name == "three_phase_fake_wells"
        inj = 1
        prod = nc
        G.grid.pore_volumes[inj] *= 1000
        G.grid.pore_volumes[prod] *= 1000

        bar = 1e5
        p0 = 100*bar # 100 bar

        mu = 1e-3    # 1 cP
        cl = 1e-5/bar
        pRef = 100*bar
        rhoLS = 1000
        A = AqueousPhase()
        L = LiquidPhase()
        V = VaporPhase()
        sys = ImmiscibleSystem([A, L, V])
        model = SimulationModel(G, sys, context = context)

        kr = BrooksCoreyRelPerm(sys, [2, 2, 2])
        s = model.secondary_variables
        s[:RelativePermeabilities] = kr
        s[:PhaseMassDensities] = ConstantCompressibilityDensities(sys, pRef, rhoLS, cl)

        tot_time = sum(timesteps)
        forces = setup_forces(model)

        p_init = repeat([p0], nc)
        p_init[inj] = 2*p0
        p_init[prod] = p0/2

        s_inj = [0.5, 0.3, 0.2]
        s_res = [0.3, 0.3, 0.4]

        s_inj = [0.1, 0.0, 0.9]
        s_res = [0.9, 0.0, 0.1]

        s_inj = [0.1, 0.9, 0.0]
        s_res = [0.9, 0.1, 0.0]

        s_init = repeat(s_inj, 1, nc)
        s_init[:, inj] .= s_inj

        # State is dict with pressure in each cell
        init = Dict(:Pressure => p_init, :Saturations => s_init)
    elseif case_name == "spe1_like_fake_wells"
        inj = 1
        prod = nc
        G.grid.pore_volumes[inj] *= 1000
        G.grid.pore_volumes[prod] *= 1000
        # Blackoil specifics
        setup = blackoil_bench_pvt(:spe1)
        pvt = setup[:pvt]
        pvto = pvt[2]
        rhoS = setup[:rhoS]
        phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
        sat_table = saturated_table(pvto)
        sys = StandardBlackOilSystem(rs_max = sat_table, phases = phases, reference_densities = rhoS)
        model = SimulationModel(G, sys, context = context)
        forces = setup_forces(model)
        # Values inside table range
        p0 = 1.38909e7
        p_init = repeat([p0], nc)
        p_init[inj] = 3.45751e7
        p_init[prod] = 1.82504e6

        sw = 0.1
        sg = zeros(nc)
        sg[inj] = 0.9
        bo = map( (P, Sg) -> BlackOilX(sys, P, sw = sw, rs = 50.0, sg = Sg), p_init, sg)
        init = Dict(:Pressure => p_init, :ImmiscibleSaturation => sw, :BlackOilUnknown => bo)
    elseif case_name == "simple_compositional_fake_wells"
        inj = 1
        prod = nc
        G.grid.pore_volumes[inj] *= 1000
        G.grid.pore_volumes[prod] *= 1000
        co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
        c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)
        c10 = MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488)
        

        z0 = [0.5, 0.3, 0.2]
        zi = [0.99, 0.01-1e-3, 1e-3]
        mixture = MultiComponentMixture([co2, c1, c10], names = ["CO2", "C1", "C10"])

        p0 = 75e5
        T0 = 423.25

        n = length(z0)
        eos = GenericCubicEOS(mixture)
        nc = number_of_cells(G)
        L, V = LiquidPhase(), VaporPhase()
        # Define system and realize on grid
        sys = MultiPhaseCompositionalSystemLV(eos, (L, V))
        model = SimulationModel(G, sys, context = context)

        kr = BrooksCoreyRelPerm(sys, [2, 3])
        s = model.secondary_variables
        s[:RelativePermeabilities] = kr
        parameters = setup_parameters(model, Temperature = T0)

        tot_time = sum(timesteps)
        forces = setup_forces(model)

        p_init = repeat([p0], nc)
        p_init[inj] = 2*p0
        p_init[prod] = p0/2

        z_init = repeat(z0, 1, nc)
        z_init[:, inj] .= zi

        # State is dict with pressure in each cell
        init = Dict(:Pressure => p_init, :OverallMoleFractions => z_init)
    elseif case_name == "compositional_three_phases"
        inj = 1
        prod = nc
        G.grid.pore_volumes[inj] *= 1000
        G.grid.pore_volumes[prod] *= 1000
        co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
        c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)
        c10 = MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488)

        z0 = [0.5, 0.3, 0.2]
        zi = [0.99, 0.01-1e-3, 1e-3]
        mixture = MultiComponentMixture([co2, c1, c10], names = ["CO2", "C1", "C10"])

        p0 = 75e5
        T0 = 423.25

        n = length(z0)
        eos = GenericCubicEOS(mixture)
        nc = number_of_cells(G)
        L, V, A = LiquidPhase(), VaporPhase(), AqueousPhase()
        # Define system and realize on grid
        sys = MultiPhaseCompositionalSystemLV(eos, (A, L, V))
        model = SimulationModel(G, sys, context = context)

        kr = BrooksCoreyRelPerm(sys, [2, 3, 2])
        s = model.secondary_variables
        s[:RelativePermeabilities] = kr
        parameters = setup_parameters(model, Temperature = T0)

        tot_time = sum(timesteps)
        forces = setup_forces(model)

        p_init = repeat([p0], nc)
        p_init[inj] = 2*p0
        p_init[prod] = p0/2

        z_init = repeat(z0, 1, nc)
        z_init[:, inj] .= zi

        sw_init = repeat([0.2], nc)
        sw_init[inj] = 0.5

        # State is dict with pressure in each cell
        init = Dict(:Pressure => p_init, :OverallMoleFractions => z_init, :ImmiscibleSaturation => sw_init)
    else
        error("Unknown case $case_name")
    end
    # Model parameters
    if isnothing(parameters)
        parameters = setup_parameters(model)
    end
    state0 = setup_state(model, init)
    return (state0, model, parameters, forces, timesteps)
end
