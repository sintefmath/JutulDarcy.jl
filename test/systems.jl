using Test, JutulDarcy, Jutul, MultiComponentFlash

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
            allocs_c = @allocations JutulDarcy.number_of_components(sys)
            allocs_ph = @allocations JutulDarcy.number_of_phases(sys)
        end
    end
end
