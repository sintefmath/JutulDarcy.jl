using Jutul, JutulDarcy, Test, LinearAlgebra
function solve_thermal(;
        nc = 10,
        time = 1000.0,
        nstep = 100,
        poro = 0.1,
        perm = 9.8692e-14,
        use_blocks = false
    )
    T = time
    tstep = repeat([T/nstep], nstep)
    G = get_1d_reservoir(nc, poro = poro, perm = perm)
    nc = number_of_cells(G)

    G[:porosity][1] *= 1000
    G[:porosity][end] *= 1000
    bar = 1e5
    p0 = repeat([1000*bar], nc)
    p0[1] = 2000*bar
    p0[end] = 500*bar

    s0 = zeros(2, nc)
    s0[2, :] .= 1.0
    s0[1, 1] = 1.0
    s0[1, 2] = 0.0

    # Define system and realize on grid
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    D = discretized_domain_tpfv_flow(G)
    if use_blocks
        l = BlockMajorLayout()
    else
        l = EquationMajorLayout()
    end
    ctx = DefaultContext(matrix_layout = l)

    model = SimulationModel(D, sys, data_domain = G, context = ctx)
    JutulDarcy.add_thermal_to_model!(model)
    push!(model.output_variables, :Temperature)
    kr = BrooksCoreyRelativePermeabilities(2, [2.0, 2.0])
    replace_variables!(model, RelativePermeabilities = kr)
    forces = setup_forces(model)

    parameters = setup_parameters(model,
                                PhaseViscosities = [1e-3, 1e-3],
                                RockThermalConductivities = 1e-2,
                                FluidThermalConductivities = 1e-2,
                                RockDensity = 1e3,
                                ComponentHeatCapacity = 10000.0,
                                RockHeatCapacity = 500.0)
    T0 = repeat([273.15], nc)
    T0[1] = 500.0
    state0 = setup_state(model, Pressure = p0, Saturations = s0, Temperature = T0)

    states, reports = simulate(state0, model, tstep,
        parameters = parameters, forces = forces, info_level = -1)
    pushfirst!(states, state0)
    return states, reports
end

using Test
@testset "simple_thermal" begin
    for use_blocks in [true, false]
        states, = solve_thermal(nc = 10, use_blocks = use_blocks);
        T = states[end][:Temperature]
        # Check that the first cell fullfills BC
        @test 490 < T[2] < 500
        # Check monotone temp profile
        @test all(x -> x < 0, diff(T))
    end
end
##
function solve_thermal_wells(;
        nx = 10,
        ny = nx,
        nz = 2,
        block_backend = false,
        thermal = true,
        simple_well = false,
        single_phase = false,
        energy_formulation = :thermal
    )
    day = 3600*24
    bar = 1e5
    g = CartesianMesh((nx, ny, nz), (1000.0, 1000.0, 100.0))
    Darcy = 9.869232667160130e-13
    K = repeat([0.1*Darcy], 1, number_of_cells(g))
    res = reservoir_domain(g, porosity = 0.1, permeability = K)
    # Vertical well in (1, 1, *), producer in (nx, ny, 1)
    P = setup_vertical_well(res, 1, 1, name = :Producer, simple_well = simple_well)
    I = setup_well(res, [(nx, ny, 1)], name = :Injector, simple_well = simple_well)
    rhoWS = 1000.0
    rhoGS = 700.0
    if single_phase
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
        nph = 1
        i_mix = [1.0]
        rhoS = [rhoWS]
        c = [1e-6/bar]
    else
        # Set up a two-phase immiscible system
        phases = (AqueousPhase(), VaporPhase())
        rhoS = [rhoWS, rhoGS]
        sys = ImmiscibleSystem(phases, reference_densities = rhoS)
        nph = 2
        i_mix = [0.0, 1.0]
        c = [1e-6/bar, 1e-5/bar]
    end
    wells = [I, P]
    model, parameters = setup_reservoir_model(res, sys,
        thermal = thermal,
        energy_formulation = energy_formulation,
        wells = wells,
        block_backend = block_backend,
        extra_out = true
    )
    # Replace the density function with our custom version for wells and reservoir
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    replace_variables!(model, PhaseMassDensities = ρ)
    dt = repeat([30.0]*day, 12*5)
    rate_target = TotalRateTarget(sum(parameters[:Reservoir][:FluidVolume])/sum(dt))
    bhp_target = BottomHolePressureTarget(50*bar)

    ictrl = InjectorControl(rate_target, i_mix, density = rhoGS, temperature = 300.0)

    pctrl = ProducerControl(bhp_target)

    controls = Dict(:Injector => ictrl,
                    :Producer => pctrl)
    forces = setup_reservoir_forces(model, control = controls)
    if thermal
        state0 = setup_reservoir_state(model,
            Pressure = 150*bar,
            Saturations = [1.0, 0.0],
            Temperature = 300.0
        )
    else
        state0 = setup_reservoir_state(model,
            Pressure = 150*bar,
            Saturations = [1.0, 0.0]
        )
    end
    result = simulate_reservoir(state0, model, dt,
        forces = forces, parameters = parameters, info_level = -1)
    return (result.states, result.result.reports, dt, model)
end

@testset "thermal wells" begin
    @testset "well types and backends" begin
        for simple_well in [true, false]
            for block_backend in [false, true]
                if block_backend
                    bs = "block"
                else
                    bs = "scalar"
                end
                if simple_well
                    ws = "simple"
                else
                    ws = "ms"
                end
                @testset "$bs $ws well" begin
                    states, reports, dt, = solve_thermal_wells(nx = 3, nz = 1,
                        thermal = true,
                        simple_well = simple_well,
                        block_backend = block_backend
                    );
                    @test length(states) == length(dt)
                end
            end
        end
    end
end

function solve_thermal_column(; energy_formulation = :total, nz = 10, height = 1000.0)
    bar = 1e5
    day = 3600*24
    g = CartesianMesh((1, 1, nz), (10.0, 10.0, height))
    res = reservoir_domain(g, porosity = 0.3, permeability = 1e-12)
    rhoS = [1000.0, 100.0]
    sys = ImmiscibleSystem((AqueousPhase(), VaporPhase()), reference_densities = rhoS)
    model, parameters = setup_reservoir_model(res, sys,
        thermal = true,
        energy_formulation = energy_formulation,
        extra_out = true
    )
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = [1e-6/bar, 1e-5/bar])
    replace_variables!(model, PhaseMassDensities = ρ)
    # Heavy phase on top of light phase: Gravity segregation in a closed column
    s = zeros(2, nz)
    s[1, 1:(nz ÷ 2)] .= 1.0
    s[2, :] .= 1.0 .- s[1, :]
    state0 = setup_reservoir_state(model,
        Pressure = 100*bar,
        Saturations = s,
        Temperature = 300.0
    )
    dt = fill(10.0*day, 50)
    forces = setup_reservoir_forces(model)
    result = simulate_reservoir(state0, model, dt,
        forces = forces, parameters = parameters, info_level = -1)
    return (result, model, parameters)
end

function solve_deep_injector(energy_formulation)
    bar = 1e5
    day = 3600*24
    nz = 20
    g = CartesianMesh((1, 1, nz), (100.0, 100.0, 100.0), origin = [0.0, 0.0, 2000.0])
    res = reservoir_domain(g, porosity = 0.2, permeability = 1e-12)
    W = setup_well(res, [(1, 1, nz)], name = :Injector, simple_well = false,
        reference_depth = 0.0, WIth = 0.0)
    sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
    model, parameters = setup_reservoir_model(res, sys,
        wells = [W],
        thermal = true,
        energy_formulation = energy_formulation,
        extra_out = true
    )
    state0 = setup_reservoir_state(model, Pressure = 200*bar, Temperature = 300.0)
    ctrl = InjectorControl(TotalRateTarget(0.01), [1.0], density = 1000.0, temperature = 300.0)
    forces = setup_reservoir_forces(model, control = Dict(:Injector => ctrl))
    result = simulate_reservoir(state0, model, fill(1.0*day, 30),
        forces = forces, parameters = parameters, info_level = -1)
    T_well = result.result.states[end][:Injector][:Temperature]
    return T_well[end] - T_well[1]
end

@testset "total energy formulation" begin
    @testset "closed system conserves total energy" begin
        result, model, = solve_thermal_column(energy_formulation = :total)
        @test JutulDarcy.model_energy_formulation(model) == :total
        rstates = result.result.states
        E(state, k) = sum(state[:Reservoir][k])
        E0_tot = E(rstates[1], :TotalEnergy)
        E_tot = E(rstates[end], :TotalEnergy)
        E0_th = E(rstates[1], :TotalThermalEnergy)
        E_th = E(rstates[end], :TotalThermalEnergy)
        # Potential energy released by gravity segregation is converted to
        # thermal energy, while the total energy is conserved (up to the
        # nonlinear solver tolerance).
        ΔE_th = E_th - E0_th
        ΔE_pot = (E_tot - E_th) - (E0_tot - E0_th)
        @test ΔE_th > 0
        @test ΔE_pot < 0
        @test abs(E_tot - E0_tot) < 0.01*abs(ΔE_pot)
        @test minimum(result.states[end][:Temperature]) > 300.0
        # Thermal formulation: Thermal energy is conserved instead, and the
        # change is much smaller than what is released as potential energy.
        result, model, = solve_thermal_column(energy_formulation = :thermal)
        @test JutulDarcy.model_energy_formulation(model) == :thermal
        rstates = result.result.states
        @test abs(E(rstates[end], :TotalThermalEnergy) - E(rstates[1], :TotalThermalEnergy)) < 0.01*abs(ΔE_pot)
    end
    @testset "deep injector" begin
        # Injection through an insulated 2 km multisegment well. Conservation of
        # thermal energy gives spurious cooling as the pressure increases
        # downwards, approximately g*Δz/c. With the total energy formulation,
        # the work done by gravity is accounted for.
        ΔT = Dict()
        for form in (:thermal, :total)
            ΔT[form] = solve_deep_injector(form)
        end
        @test ΔT[:thermal] < -3.0
        @test abs(ΔT[:total]) < 1.0
    end
    @testset "wells" begin
        for simple_well in [true, false]
            @testset "simple=$simple_well" begin
                states, reports, dt, model = solve_thermal_wells(nx = 3, nz = 2,
                    simple_well = simple_well,
                    energy_formulation = :total
                )
                @test length(states) == length(dt)
                for k in (:Reservoir, :Injector, :Producer)
                    @test JutulDarcy.model_energy_formulation(model[k]) == :total
                end
                @test haskey(states[end], :TotalEnergy)
                @test model.models[:Producer].equations[:energy_conservation] isa Jutul.ConservationLaw{:TotalEnergy}
            end
        end
        @test_throws ArgumentError solve_thermal_wells(nx = 3, nz = 1,
            energy_formulation = :kinetic
        )
    end
end

function solve_isothermal_injection(sys, T0; mix = [1.0], density = 1000.0, state_arg = NamedTuple(), kwarg...)
    bar = 1e5
    day = 3600*24
    g = CartesianMesh((5, 1, 1), (500.0, 100.0, 10.0), origin = [0.0, 0.0, 1000.0])
    d = reservoir_domain(g, permeability = 1e-13, porosity = 0.1)
    # Separate top node so that it is not directly connected to the reservoir
    I = setup_well(d, [(1, 1, 1)], name = :Injector, simple_well = false, use_top_node = true)
    P = setup_well(d, [(5, 1, 1)], name = :Producer, simple_well = false)
    model, parameters = setup_reservoir_model(d, sys, wells = [I, P], extra_out = true; kwarg...)
    state0 = setup_reservoir_state(model; Pressure = 150bar, Temperature = T0, state_arg...)
    ctrl = Dict(
        :Injector => InjectorControl(TotalRateTarget(1e-3), mix, density = density, temperature = T0),
        :Producer => ProducerControl(BottomHolePressureTarget(140bar))
    )
    forces = setup_reservoir_forces(model, control = ctrl)
    result = simulate_reservoir(state0, model, fill(10.0*day, 20),
        forces = forces, parameters = parameters, info_level = -1)
    return (result, model)
end

@testset "tabulated internal energy" begin
    tables = JutulDarcy.Geothermal.geothermal_setup_tables()
    for T0 in [300.0, 350.0]
        result, model = solve_isothermal_injection(:geothermal, T0)
        for k in (:Reservoir, :Injector, :Producer)
            m = model.models[k]
            @test Jutul.get_secondary_variables(m)[:FluidInternalEnergy] isa JutulDarcy.PressureTemperatureDependentInternalEnergy
            # Heat capacity is not used and should not be present
            @test !haskey(Jutul.get_secondary_variables(m), :ComponentHeatCapacity)
            @test !haskey(Jutul.get_parameters(m), :ComponentHeatCapacity)
        end
        states = result.result.states
        # Internal energy is taken from the table
        res = states[end][:Reservoir]
        p, T = res[:Pressure][1], res[:Temperature][1]
        @test res[:FluidInternalEnergy][1] ≈ tables[:internal_energy](p, T) rtol = 1e-8
        # Injected enthalpy is evaluated from the same tables at the injection
        # temperature, so the top node of the injector is at that temperature.
        inj = states[end][:Injector]
        p_top = inj[:Pressure][1]
        H_inj = tables[:internal_energy](p_top, T0) + p_top/tables[:density](p_top, T0)
        @test inj[:FluidEnthalpy][1] ≈ H_inj rtol = 1e-8
        @test inj[:Temperature][1] ≈ T0 atol = 1e-4
        # Injection at reservoir temperature: Only small changes due to
        # Joule-Thomson effects.
        @test all(x -> abs(x - T0) < 0.5, res[:Temperature])
    end
    @testset "table edges" begin
        U = JutulDarcy.PressureTemperatureDependentInternalEnergy(tables[:internal_energy])
        T_lo, T_hi = first(tables[:internal_energy].Y), last(tables[:internal_energy].Y)
        u(T) = JutulDarcy.tabulated_internal_energy(U, 5e5, T, 1)
        dudT(T) = (u(T + 0.01) - u(T - 0.01))/0.02
        # The tables are padded with a constant first column. This should not
        # give zero heat capacity inside the table, and the internal energy
        # should be extrapolated linearly outside it.
        for T in [T_lo - 20.0, T_lo + 0.5, 273.15, T_hi + 20.0]
            @test dudT(T) > 1000.0
        end
        # Temperature limit is lowered so that sub-zero iterates are allowed
        model = setup_reservoir_model(reservoir_domain(CartesianMesh((2, 1, 1))), :geothermal)
        @test Jutul.get_primary_variables(model[:Reservoir])[:Temperature].min == 200.0
    end
    @testset "co2-brine $physics" for physics in (:kvalue, :immiscible)
        T0 = 320.0
        result, model = solve_isothermal_injection(:co2brine, T0,
            mix = [0.0, 1.0],
            density = 1.87,
            thermal = true,
            co2_physics = physics,
            state_arg = physics == :kvalue ? (OverallMoleFractions = [1.0, 0.0], ) : (Saturations = [1.0, 0.0], )
        )
        rmodel = model.models[:Reservoir]
        @test Jutul.get_secondary_variables(rmodel)[:FluidInternalEnergy] isa JutulDarcy.PressureTemperatureDependentInternalEnergy
        @test !haskey(Jutul.get_secondary_variables(rmodel), :ComponentHeatCapacity)
        T = result.result.states[end][:Reservoir][:Temperature]
        @test all(x -> abs(x - T0) < 0.5, T)
    end
end
