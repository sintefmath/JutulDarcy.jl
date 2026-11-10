using Test, JutulDarcy, Jutul, MultiComponentFlash

@testset "Black-oil water-only startup" begin
    deck = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE9", "SPE9.DATA")
    props = JutulDarcy.GeoEnergyIO.parse_data_file(deck)["PROPS"]
    pvtw = JutulDarcy.deck_pvt_water(props)
    pvto = JutulDarcy.deck_pvt_oil(props)
    pvdg = JutulDarcy.deck_pvt_gas(props)
    oil, water, gas = only(props["DENSITY"])
    darcy, bar, day = si_units(:darcy, :bar, :day)
    domain = reservoir_domain(reservoir_mesh((1, 1, 16), (1.0, 1.0, 16.0));
        permeability = darcy, porosity = 0.2)
    model, parameters = JutulDarcy.setup_reservoir_model_from_blackoil_tables(domain;
        pvtw, pvto, pvdg, reference_densities = [water, oil, gas], extra_out = true,
        extra_outputs = [:Saturations, :Rs, :TotalMasses])
    reservoir = reservoir_model(model)
    kr = BrooksCoreyRelativePermeabilities(reservoir.system, 2.0)
    set_secondary_variables!(reservoir; RelativePermeabilities = kr)
    sw = [ones(8); zeros(8)]
    saturations = vcat(sw', (1 .- sw)', zeros(1, length(sw)))
    state0 = setup_reservoir_state(model; Pressure = fill(100*bar, 16), Saturations = saturations, Rs = zeros(16))
    unknowns = state0[:Reservoir][:BlackOilUnknown]
    @test all(x -> x.phases_present == JutulDarcy.OilAndGas && iszero(x.val), unknowns[1:8])
    @test all(x -> x.phases_present == JutulDarcy.OilOnly && iszero(x.val), unknowns[9:16])
    initial = Jutul.evaluate_all_secondary_variables(reservoir, state0[:Reservoir], parameters[:Reservoir])
    @test initial[:Saturations] ≈ saturations
    @test all(iszero, initial[:TotalMasses][2:3, 1:8])
    case = JutulCase(model, [0.05*day, 0.5*day]; state0, parameters)
    result = simulate_reservoir(case; info_level = -1, max_timestep_cuts = 0)
    @test length(result.states) == length(case.dt)
    final = last(result.states)
    @test final[:Saturations][1, 8] < 1
    @test sum(final[:TotalMasses]; dims = 2) ≈ sum(initial[:TotalMasses]; dims = 2) rtol = 1e-6
    @test maximum(abs, final[:Rs]) < 1e-6
end

@testset "Black-oil water-filled saturation coordinates" begin
    # Water-filled cells use gas saturation as their coordinate, including
    # when a small, nonzero hydrocarbon inventory was specified.
    tab = Jutul.get_1d_interpolator([1e6, 2e7], [10.0, 100.0])
    eps_s = JutulDarcy.MINIMUM_COMPOSITIONAL_SATURATION
    # Stay strictly inside the cutoff; equality uses normal phase switching.
    sw = 1 - 5*eps_s
    rem = 1 - sw
    for (rs_tab, rv_tab) in ((tab, nothing), (nothing, tab), (tab, tab))
        for sg in (0.0, rem/4, rem)
            x = JutulDarcy.blackoil_unknown_init(rs_tab, rv_tab, sw, rem-sg, sg, 0.0, 0.0, 1e7)
            @test x.phases_present == JutulDarcy.OilAndGas
            @test x.val ≈ sg
        end
    end
    for sg in (0.0, rem/4, rem), dx in (-0.1, 0.0, 0.1)
        x = [JutulDarcy.BlackOilX(sg)]
        JutulDarcy.varswitch_update_inner!(x, 1, dx, Inf, 0.2, tab, nothing,
            false, true, [1e7], sw, (1e-10, 1e-10, eps_s), 1.0)
        @test x[1].phases_present == JutulDarcy.OilAndGas
        @test x[1].val ≈ clamp(sg + dx, 0.0, rem)
    end
    @testset "Cutoff boundary" begin
        for eps_s in (1e-12, 1e-10, 1e-8)
            cutoff = 1 - 10*eps_s
            @test !JutulDarcy.blackoil_water_filled(prevfloat(cutoff), eps_s)
            @test !JutulDarcy.blackoil_water_filled(cutoff, eps_s)
            @test JutulDarcy.blackoil_water_filled(nextfloat(cutoff), eps_s)
        end
    end
    @testset "Initialization uses the system tolerance" begin
        for eps_s in (1e-12, 1e-10, 1e-8), (rs_tab, rv_tab) in ((tab, nothing), (nothing, tab), (tab, tab))
            sys = StandardBlackOilSystem(; rs_max = rs_tab, rv_max = rv_tab, eps_s)
            for (factor, water_filled) in ((5, true), (20, false))
                sw = 1 - factor*eps_s
                rem = 1 - sw
                sg = isnothing(rs_tab) ? rem : 0.0
                so = rem - sg
                single_phase = isnothing(rs_tab) ? JutulDarcy.GasOnly : JutulDarcy.OilOnly
                expected_phase = water_filled ? JutulDarcy.OilAndGas : single_phase
                expected_value = water_filled ? sg : 0.0
                direct = JutulDarcy.blackoil_unknown_init(rs_tab, rv_tab, sw, so, sg, 0.0, 0.0, 1e7; eps_s)
                initialized = BlackOilX(sys, 1e7; sw, so, sg, rs = 0.0, rv = 0.0)
                for x in (direct, initialized)
                    @test x.phases_present == expected_phase
                    @test x.val ≈ expected_value
                end
            end
        end
    end
end

@testset "Black-oil optimization initial state" begin
    rs_tab = Jutul.get_1d_interpolator([1e6, 2e7], [10.0, 100.0])
    rv_tab = Jutul.get_1d_interpolator([1e6, 2e7], [0.001, 0.01])
    for water in (false, true), eps_s in (water ? (1e-12, 1e-10, 1e-8) : (1e-10,)), (rs, rv) in ((rs_tab, nothing), (nothing, rv_tab), (rs_tab, rv_tab))
        phases = water ? (AqueousPhase(), LiquidPhase(), VaporPhase()) : (LiquidPhase(), VaporPhase())
        sys = StandardBlackOilSystem(; rs_max = rs, rv_max = rv, phases, eps_s)
        domain = reservoir_domain(reservoir_mesh((3, 1), (3.0, 1.0));
            permeability = si_unit(:darcy), porosity = 0.2)
        model, parameters = setup_reservoir_model(domain, sys; extra_out = true)
        sw = water ? [1.0, 1-5*sys.s_eps, 0.2] : zeros(3)
        sg = water ? [0.0, isnothing(rs) ? 1-sw[2] : 0.0, 0.3] : [0.0, 0.3, 1.0]
        so = 1 .- sw .- sg
        sat = water ? vcat(sw', so', sg') : vcat(so', sg')
        ratios = Dict{Symbol, Any}()
        isnothing(rs) || (ratios[:Rs] = zeros(3))
        isnothing(rv) || (ratios[:Rv] = zeros(3))
        state0 = setup_reservoir_state(model; Pressure = fill(1e7, 3), Saturations = sat, ratios...)
        if water
            @test state0[:Reservoir][:BlackOilUnknown][2].phases_present == JutulDarcy.OilAndGas
        end
        case = JutulCase(model, [1.0]; state0, parameters)
        dprm = setup_reservoir_dict_optimization(case; verbose = false)
        extracted = dprm.parameters[:state0]
        @test extracted[:LiquidSaturation] ≈ so
        @test extracted[:VaporSaturation] ≈ sg
        if !isnothing(rs)
            expected = water ? [0.0, 0.0, rs(1e7)] : [0.0, rs(1e7), rs(1e7)]
            @test extracted[:Rs] ≈ expected
        end
        if !isnothing(rv)
            expected = water ? [0.0, 0.0, rv(1e7)] : [rv(1e7), rv(1e7), 0.0]
            @test extracted[:Rv] ≈ expected
        end
        rebuilt = dprm.setup_function(dprm.parameters)
        rmodel = reservoir_model(rebuilt.model)
        restored = similar(sat)
        Jutul.update_secondary_variable!(restored, rmodel[:Saturations], rmodel,
            Jutul.JutulStorage(rebuilt.state0[:Reservoir]), axes(sat, 2))
        @test restored ≈ sat
    end
end

@testset "number of components and phases" begin
    @testset "immiscible" begin
        for phases in [
                (LiquidPhase(),),
                (LiquidPhase(), VaporPhase()),
                (LiquidPhase(), AqueousPhase(), VaporPhase())
            ]
            sys = ImmiscibleSystem(phases)
            nc = JutulDarcy.number_of_components(sys)
            nph = JutulDarcy.number_of_phases(sys)

            @test nph == length(phases) == nc
            @test typeof(sys).parameters[end] ==
                JutulDarcy.get_reference_phase_index(sys)
            allocs_c = @allocations JutulDarcy.number_of_components(sys)
            allocs_ph = @allocations JutulDarcy.number_of_phases(sys)

            @test allocs_c == 0
            @test allocs_ph == 0
        end
    end

    @testset "compositional" begin
        for include_water in [true, false]
            setup = JutulDarcy.blackoil_bench_pvt(:spe1)
            pvt = setup[:pvt]
            pvto = pvt[2].tab[1]
            rhoS = setup[:rhoS]

            if include_water
                phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
            else
                phases = (LiquidPhase(), VaporPhase())
            end
            sat_table = get_1d_interpolator(pvto.sat_pressure, pvto.rs, cap_end = false)
            sys = StandardBlackOilSystem(rs_max = sat_table, phases = phases, reference_densities = rhoS)

            nc = JutulDarcy.number_of_components(sys)
            nph = JutulDarcy.number_of_phases(sys)

            @test nph == length(phases) == nc
            @test typeof(sys).parameters[end] ==
                JutulDarcy.get_reference_phase_index(sys)
            allocs_c = @allocations JutulDarcy.number_of_components(sys)
            allocs_ph = @allocations JutulDarcy.number_of_phases(sys)
        end
    end

    @testset "compositional" begin
        for include_water in [true, false]

            if include_water
                phases = [AqueousPhase(), LiquidPhase(), VaporPhase()]
            else
                phases = [LiquidPhase(), VaporPhase()]
            end
            co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
            c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)
            c10 = MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488)
            components = [co2, c1, c10]

            mixture = MultiComponentMixture(components, names = ["CO2", "C1", "C10"])
            eos = GenericCubicEOS(mixture)
            sys = MultiPhaseCompositionalSystemLV(eos, phases)

            nc = JutulDarcy.number_of_components(sys)
            nph = JutulDarcy.number_of_phases(sys)

            @test nph == length(phases)
            @test nc == length(components) + include_water
            @test typeof(sys).parameters[end] ==
                JutulDarcy.get_reference_phase_index(sys)
            overridden = MultiPhaseCompositionalSystemLV(eos, phases;
                reference_phase_index = 1)
            @test JutulDarcy.get_reference_phase_index(overridden) == 1
            @test typeof(overridden).parameters[end] == 1
            allocs_c = @allocations JutulDarcy.number_of_components(sys)
            allocs_ph = @allocations JutulDarcy.number_of_phases(sys)
        end
    end
end

@testset "Reference phase overrides" begin
    single_phase = SinglePhaseSystem()
    @test JutulDarcy.get_reference_phase_index(single_phase) == 1
    @test typeof(single_phase).parameters[end] == 1
    phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
    immiscible = ImmiscibleSystem(phases; reference_phase_index = 1)
    black_oil = StandardBlackOilSystem(
        phases = phases, reference_phase_index = 1)
    well = JutulDarcy.SimpleWellSystem(immiscible)
    for system in (immiscible, black_oil, well)
        @test JutulDarcy.get_reference_phase_index(system) == 1
        @test typeof(system).parameters[end] == 1
    end
end
