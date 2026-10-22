using Jutul, JutulDarcy, MultiComponentFlash, Test

function compositional_keyword_deck(; water = false, equil = false)
    water_props = water ? "DENSITY\n800 1000 1 /\nPVTW\n75 1 0 1 0 /\nSWOF\n0 0 1 0\n1 1 0 0 /" : ""
    solution = equil ? "EQUIL\n1010 75 1500 0 900 0 2* 0 1 /" : "PRESSURE\n2*75 /\nSGAS\n2*0.298017 /\nZMF\n2*0.6 2*0.1 2*0.3 /\nTEMPI\n2*150 /"
    return """
RUNSPEC
METRIC
OIL
GAS
$(water ? "WATER" : "")
COMPS
3 /
DIMENS
2 1 1 /
GRID
DXV
2*100 /
DYV
100 /
DZV
20 /
TOPS
2*1000 /
PORO
2*0.25 /
PERMX
2*100 /
PERMY
2*100 /
PERMZ
2*10 /
PROPS
EOS
PR /
CNAMES
DECANE CO2 METHANE /
MW
142.28 44 16.04 /
PCRIT
21.03 73.773 45.992 /
TCRIT
617.7 304.128 190.564 /
VCRIT
0.6098 0.09412 0.09863 /
ACF
0.4884 0.22394 0.01142 /
OMEGAA
0.46 0.45 0.47 /
OMEGAB
0.078 0.076 0.079 /
SSHIFT
0.01 0.02 0.03 /
SSHIFTS
0.04 0.05 0.06 /
LBCCOEF
0.11 0.024 0.059 -0.041 0.0094 /
RTEMP
150 /
COMPVD
1000 0.6 0.1 0.3 1 70 /
$water_props
SGOF
0 0 1 0
1 1 0 0 /
SOLUTION
$solution
$(water && !equil ? "SWAT\n2*0 /" : "")
SCHEDULE
WELSPECS
I CHILD 1 1 1010 GAS /
/
GRUPTREE
CHILD PARENT /
PARENT FIELD /
/
COMPDAT
I 1 1 1 1 OPEN 2* 0.2 /
/
WELLSTRE
A 0 0.9 0.1 /
B 0 0.1 0.9 /
/
GINJGAS
FIELD STREAM A /
/
WCONINJE
I GAS OPEN RATE 10 1* 100 /
/
TSTEP
0.001 /
GINJGAS
PARENT STREAM B /
/
TSTEP
0.001 /
END
"""
end

@testset "Compositional deck keyword physics" begin
    MCF = MultiComponentFlash
    mktempdir() do dir
        path = joinpath(dir, "COMPOSITIONAL.DATA")
        write(path, compositional_keyword_deck())
        data = JutulDarcy.GeoEnergyIO.parse_data_file(path; silent = true)
        sys, = JutulDarcy.parse_physics_types(data)
        eos = sys.equation_of_state
        facility_sys, = JutulDarcy.parse_physics_types(data; for_facility = true)
        facility_eos = facility_sys.equation_of_state
        base = GenericCubicEOS(eos.mixture, PengRobinson())
        cond = (p = 7.5e6, T = 423.15, z = [0.6, 0.1, 0.3])
        for i in 1:3
            @test MCF.A_i(eos, cond, i)/MCF.A_i(base, cond, i) ≈ data["PROPS"]["OMEGAA"][1][i]/base.ω_a
            @test MCF.B_i(eos, cond, i)/MCF.B_i(base, cond, i) ≈ data["PROPS"]["OMEGAB"][1][i]/base.ω_b
        end
        @test eos.volume_shift == (0.01, 0.02, 0.03)
        @test facility_eos.volume_shift == (0.04, 0.05, 0.06)
        for keyword in ("OMEGAA", "OMEGAB")
            props = Dict(keyword => data["PROPS"][keyword])
            parameters = JutulDarcy.deck_cubic_coefficients(PengRobinson(), props, 3)
            if keyword != "OMEGAA"
                @test parameters.omega_a == (base.ω_a, base.ω_a, base.ω_a)
            end
            if keyword != "OMEGAB"
                @test parameters.omega_b == (base.ω_b, base.ω_b, base.ω_b)
            end
        end
        for F in (Float64, Float32)
            static = MCF.make_eos_immutable(eos; float_type = F)
            @test isbitstype(typeof(static))
            @test eltype(static.type.omega_a) == F
            @test MCF.force_coefficients(static, cond).A_i ≈ MCF.force_coefficients(eos, cond).A_i rtol = 1e-6
        end
        flash = MCF.flashed_mixture_2ph(eos, cond)
        @test MCF.molar_volume(eos, cond.p, cond.T, flash.liquid) != MCF.molar_volume(facility_eos, cond.p, cond.T, flash.liquid)
        coefficients = Tuple(data["PROPS"]["LBCCOEF"])
        mu = JutulDarcy.compositional_lbc_viscosities(eos, cond.p, cond.T, flash, coefficients)
        @test mu.μ_l ≈ MCF.lbc_viscosity(eos, cond.p, cond.T, flash.liquid; coeff = coefficients)
        @test mu.μ_l != MCF.lbc_viscosity(eos, cond.p, cond.T, flash.liquid)

        streams = JutulDarcy.setup_well_streams(data["SCHEDULE"]["WELSPECS"])
        steps = data["SCHEDULE"]["STEPS"]
        JutulDarcy.parse_well_streams_for_step!(streams, steps[1], data["PROPS"])
        @test streams.wells["I"] == "A"
        JutulDarcy.parse_well_streams_for_step!(streams, steps[2], data["PROPS"])
        @test streams.wells["I"] == "B"
        override = Dict("WINJGAS" => [["I", "STREAM", "A", "Default", 0]])
        JutulDarcy.parse_well_streams_for_step!(streams, override, data["PROPS"])
        @test streams.wells["I"] == "A"
        inherit = Dict("WINJGAS" => [["I", "GRUP", "Default", "Default", 0]])
        JutulDarcy.parse_well_streams_for_step!(streams, inherit, data["PROPS"])
        @test streams.wells["I"] == "B"
        wildcard = Dict("GINJGAS" => [["PAR*", "ST", "A", "Default", 0]])
        JutulDarcy.parse_well_streams_for_step!(streams, wildcard, data["PROPS"])
        @test streams.wells["I"] == "A"
        inherited_group = Dict("GINJGAS" => [["PARENT", "GRUP", "Default", "Default", 0]])
        JutulDarcy.parse_well_streams_for_step!(streams, inherited_group, data["PROPS"])
        @test streams.wells["I"] == "A"
        unsupported = Dict("GINJGAS" => [["FIELD", "GV", "FIELD", "Default", 0]])
        @test_throws ArgumentError JutulDarcy.parse_well_streams_for_step!(streams, unsupported, data["PROPS"])

        for water in (false, true)
            write(path, compositional_keyword_deck(; water = water))
            case = setup_case_from_data_file(path)
            model = reservoir_model(case)
            well_system = case.model[:I].system
            facility_system = case.model[:Facility].system.multiphase
            @test model.system.equation_of_state.volume_shift == (0.01, 0.02, 0.03)
            @test well_system.equation_of_state.volume_shift == (0.01, 0.02, 0.03)
            @test facility_system.equation_of_state.volume_shift == (0.04, 0.05, 0.06)
            @test model.secondary_variables[:PhaseViscosities].coefficients == coefficients
            @test case.model[:I].secondary_variables[:PhaseViscosities].coefficients == coefficients
            controls = [f[:Facility].control[:I] for f in case.forces]
            @test controls[1].injection_mixture != controls[2].injection_mixture
            expected = [0.0, 0.1*0.044, 0.9*0.01604]
            expected /= sum(expected)
            @test controls[2].injection_mixture[1:3] ≈ expected atol = 1e-8
            surface_density, = JutulDarcy.select_injector_mixture_spec(facility_system, "I", streams, "GAS")
            well_density, = JutulDarcy.select_injector_mixture_spec(well_system, "I", streams, "GAS")
            @test controls[1].mixture_density ≈ surface_density
            @test surface_density != well_density
            adapted = JutulDarcy.Adapt.adapt(Jutul.DefaultContext(), well_system)
            @test collect(adapted.equation_of_state.volume_shift) ≈ [0.01, 0.02, 0.03]
            adapted_facility = JutulDarcy.Adapt.adapt(Jutul.DefaultContext(), case.model[:Facility].system)
            @test collect(adapted_facility.multiphase.equation_of_state.volume_shift) ≈ [0.04, 0.05, 0.06]
            reference = simulate_reservoir(case; info_level = -1, linear_solver = nothing)
            kernel = simulate_reservoir(case; mode = :ka, info_level = -1, linear_solver = nothing)
            @test length(reference.states) == length(kernel.states) == 2
            for (state, ref) in zip(kernel.states, reference.states)
                @test all(isfinite, state[:Pressure])
                @test state[:Pressure] ≈ ref[:Pressure] rtol = 1e-7
                @test state[:OverallMoleFractions] ≈ ref[:OverallMoleFractions] rtol = 1e-7
            end
        end
        write(path, compositional_keyword_deck(; equil = true))
        case = setup_case_from_data_file(path; skip_wells = true)
        @test case.state0[:Reservoir][:OverallMoleFractions] ≈ repeat([0.6, 0.1, 0.3], 1, 2)
        @test all(isfinite, case.state0[:Reservoir][:Pressure])

        # Surface shifts alone use the ordinary cubic EOS, with no wrapper.
        shift_only = deepcopy(data)
        delete!(shift_only["PROPS"], "OMEGAA")
        delete!(shift_only["PROPS"], "OMEGAB")
        delete!(shift_only["PROPS"], "SSHIFT")
        shift_only_reservoir, = JutulDarcy.parse_physics_types(shift_only)
        shift_only_facility, = JutulDarcy.parse_physics_types(shift_only; for_facility = true)
        @test shift_only_reservoir.equation_of_state.type isa PengRobinson
        @test shift_only_reservoir.equation_of_state.volume_shift === nothing
        @test shift_only_facility.equation_of_state.type isa PengRobinson
        @test shift_only_facility.equation_of_state.volume_shift == (0.04, 0.05, 0.06)

        delete!(data["PROPS"], "SSHIFTS")
        unchanged, = JutulDarcy.parse_physics_types(data)
        inherited, = JutulDarcy.parse_physics_types(data; for_facility = true)
        @test inherited.equation_of_state.volume_shift == unchanged.equation_of_state.volume_shift
        inherited_case = setup_case_from_data_file(data)
        @test inherited_case.model[:I].system.equation_of_state.volume_shift == (0.01, 0.02, 0.03)
        @test inherited_case.model[:Facility].system.multiphase.equation_of_state.volume_shift == (0.01, 0.02, 0.03)
        data["PROPS"]["PRCORR"] = true
        corrected, = JutulDarcy.parse_physics_types(data)
        @test corrected.equation_of_state.type.type isa PengRobinsonCorrected
    end
end

@testset "COMPVD phase interpolation" begin
    table = [1000 0.1 0.9 0 1e7; 1100 0.3 0.7 0 1e7;
             1200 0.8 0.2 1 1e7; 1300 0.6 0.4 1 1e7]
    composition = JutulDarcy.compvd_composition_function(table, 1150)
    @test composition(1050) ≈ [0.2, 0.8]
    @test composition(1125) ≈ [0.3, 0.7]
    @test composition(1175) ≈ [0.8, 0.2]
    @test composition(1250) ≈ [0.7, 0.3]
    @test composition(500) ≈ [0.1, 0.9]
    @test composition(2000) ≈ [0.6, 0.4]
    @test composition(1200, VaporPhase()) ≈ [0.3, 0.7]
    @test composition(1000, LiquidPhase()) ≈ [0.8, 0.2]
    invalid = copy(table)
    invalid[1, end-1] = 2
    @test_throws ArgumentError JutulDarcy.compvd_composition_function(invalid, 1150)
end
