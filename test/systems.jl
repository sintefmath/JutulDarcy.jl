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
    @test all(x -> x.phases_present == JutulDarcy.OilOnly && iszero(x.val),
        state0[:Reservoir][:BlackOilUnknown])
    initial = Jutul.evaluate_all_secondary_variables(reservoir, state0[:Reservoir], parameters[:Reservoir])
    case = JutulCase(model, [0.05*day, 0.5*day]; state0, parameters)
    result = simulate_reservoir(case; info_level = -1, max_timestep_cuts = 0)
    @test length(result.states) == length(case.dt)
    final = last(result.states)
    @test final[:Saturations][1, 8] < 1
    @test sum(final[:TotalMasses]; dims = 2) ≈ sum(initial[:TotalMasses]; dims = 2) rtol = 1e-6
    @test maximum(abs, final[:Rs]) < 1e-6
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
