using Jutul, JutulDarcy, GeoEnergyIO, Test

function deck_group_test_data(schedule)
    deck = """
RUNSPEC
METRIC
OIL
WATER
DIMENS
6 4 1 /
GRID
DX
24*100 /
DY
24*100 /
DZ
24*30 /
TOPS
24*1000 /
PORO
24*0.2 /
PERMX
24*200 /
PERMY
24*200 /
PERMZ
24*200 /
PROPS
DENSITY
800 1000 1 /
PVTW
200 1 1e-5 1 0 /
PVDO
100 1 1
300 1 1 /
SWOF
0 0 1 0
1 1 0 0 /
SOLUTION
PRESSURE
24*200 /
SWAT
24*0.4 /
SCHEDULE
GRUPTREE
CHILD PROD /
PROD FIELD /
INJ FIELD /
/
WELSPECS
P1 CHILD 6 1 1015 OIL /
P2 CHILD 6 4 1015 OIL /
I1 INJ 1 1 1015 WATER /
I2 INJ 1 4 1015 WATER /
/
COMPDAT
P1 6 1 1 1 OPEN 1* 1* 0.2 /
P2 6 4 1 1 OPEN 1* 1* 0.2 /
I1 1 1 1 1 OPEN 1* 1* 0.2 /
I2 1 4 1 1 OPEN 1* 1* 0.2 /
/
$schedule
END
"""
    return mktempdir() do dir
        path = joinpath(dir, "GROUP.DATA")
        write(path, deck)
        parse_data_file(path; silent = true)
    end
end

@testset "Deck group controls" begin
    schedule = """
GCONPROD
PROD ORAT 86.4 3* RATE /
/
GCONINJE
INJ WATER RATE 86.4 /
/
WCONPROD
P1 SHUT GRUP 4* 1* 100 /
P2 OPEN GRUP 4* 1* 100 /
/
WELOPEN
P1 OPEN /
/
WCONINJE
I1 WATER OPEN GRUP 2* 300 /
I2 WATER OPEN GRUP 2* 300 /
/
TSTEP
0.1 /
GCONPROD
PROD ORAT 43.2 3* RATE /
/
GCONINJE
INJ WATER RATE 43.2 /
/
TSTEP
0.1 /
WELOPEN
P1 SHUT /
I1 SHUT /
/
TSTEP
0.1 /
WELOPEN
P1 OPEN /
I1 OPEN /
/
TSTEP
0.1 /
"""
    data = deck_group_test_data(schedule)
    case = setup_case_from_data_file(data; group_use_potentials = Dict(:PROD_PRODUCTION => false, :INJ_WATER => false))
    @test Set(case.model[:Facility].domain.group_to_wells[:PROD_PRODUCTION]) == Set([:P1, :P2])
    @test case.model[:Facility].domain.lag_potentials
    @test Set(case.model[:Facility].domain.group_to_wells[:INJ_WATER]) == Set([:I1, :I2])
    @test case.forces[1][:Facility].control[:P1].target == GroupSurfaceOilRateTarget(-1e-3)
    @test case.forces[2][:Facility].control[:P1].target == GroupSurfaceOilRateTarget(-0.5e-3)
    @test case.forces[2][:Facility].control[:I1].target == GroupSurfaceWaterRateTarget(0.5e-3)
    @test isnothing(case.forces[1][:P1].mask) || any(>(0), case.forces[1][:P1].mask.values)
    @test case.forces[3][:Facility].control[:P1] isa DisabledControl
    @test case.forces[4][:Facility].control[:P1] isa ProducerControl
    result = simulate_reservoir(case; info_level = -1, timesteps = :none, error_on_incomplete = true)
    @test length(result.result.states) == 4
    for (step, state) in enumerate(result.result.states)
        f = state[:Facility]
        names = case.model[:Facility].domain.well_symbols
        p = [findfirst(isequal(w), names) for w in (:P1, :P2)]
        i = [findfirst(isequal(w), names) for w in (:I1, :I2)]
        rate = step == 1 ? 1e-3 : 0.5e-3
        @test sum(f[:SurfacePhaseRates][2, p]) ≈ -rate rtol = 1e-3
        @test sum(f[:SurfacePhaseRates][1, i]) ≈ rate rtol = 1e-3
        if step == 3
            @test all(iszero, f[:SurfacePhaseRates][:, [p[1], i[1]]])
        else
            @test f[:SurfacePhaseRates][2, p] ≈ fill(-rate/2, 2) rtol = 1e-3
            @test f[:SurfacePhaseRates][1, i] ≈ fill(rate/2, 2) rtol = 1e-3
        end
    end

    @testset "Individual requests retain group limits" begin
        # Supply finite individual oil requests while the group cap changes.
        ordinary = replace(schedule,
            "P1 SHUT GRUP 4* 1* 100" => "P1 SHUT ORAT 864 3* 1* 100",
            "P2 OPEN GRUP 4* 1* 100" => "P2 OPEN ORAT 864 3* 1* 100")
        other_case = setup_case_from_data_file(deck_group_test_data(ordinary))
        @test other_case.forces[1][:Facility].control[:P1].target isa SurfaceOilRateTarget
        @test other_case.forces[1][:Facility].limits[:P1].group_orat == -1e-3
        @test other_case.forces[2][:Facility].limits[:P1].group_orat == -0.5e-3
        other = simulate_reservoir(other_case[1:2]; info_level = -1, timesteps = :none, error_on_incomplete = true)
        names = other_case.model[:Facility].domain.well_symbols
        p = [findfirst(isequal(w), names) for w in (:P1, :P2)]
        for (step, state) in enumerate(other.result.states)
            @test sum(state[:Facility][:SurfacePhaseRates][2, p]) ≈ -(step == 1 ? 1e-3 : 0.5e-3) rtol = 1e-3
        end
    end

    @testset "WGRUPCON opt-out" begin
        excluded = replace(schedule, "P2 OPEN GRUP 4* 1* 100" => "P2 OPEN ORAT 17.28 3* 1* 100")
        excluded = "WGRUPCON\nP2 NO /\n/\n"*excluded
        excluded_case = setup_case_from_data_file(deck_group_test_data(excluded))
        @test excluded_case.model[:Facility].domain.group_to_wells[:PROD_PRODUCTION] == [:P1]
        @test !haskey(excluded_case.forces[1][:Facility].limits[:P2], :group_orat)
        other = simulate_reservoir(excluded_case[1:2]; info_level = -1, timesteps = :none, error_on_incomplete = true)
        names = excluded_case.model[:Facility].domain.well_symbols
        p1 = findfirst(isequal(:P1), names)
        p2 = findfirst(isequal(:P2), names)
        for (step, state) in enumerate(other.result.states)
            @test state[:Facility][:SurfacePhaseRates][2, p1] ≈ -(step == 1 ? 1e-3 : 0.5e-3) rtol = 1e-3
            @test state[:Facility][:SurfacePhaseRates][2, p2] ≈ -0.2e-3 rtol = 1e-3
        end
    end

    @testset "Unsupported group features" begin
        for (from, to, message) in (
                ("INJ WATER RATE 86.4", "INJ WATER VREP 3* 1", "VREP group control"),
                ("PROD ORAT 86.4 3* RATE", "PROD ORAT 86.4 3* WELL", "exceed procedures"),
                ("PROD ORAT 86.4 3* RATE", "'P*' ORAT 86.4", "wildcard group names"),
                ("PROD ORAT 86.4 3* RATE", "PROD FLD 86.4", "subordinate group constraints"),
                ("CHILD PROD /", "CHILD PROD /\nPROD CHILD /", "Cycle in GRUPTREE"))
            modified = deepcopy(data)
            # Reparse text so the parser and importer are tested together.
            source = replace(schedule, from => to)
            if from == "CHILD PROD /"
                modified["SCHEDULE"]["STEPS"][1]["GRUPTREE"] = [["CHILD", "PROD"], ["PROD", "CHILD"]]
            else
                modified = deck_group_test_data(source)
            end
            @test_throws message setup_case_from_data_file(modified)
        end
        no_group = deck_group_test_data(replace(schedule,
            "PROD ORAT 86.4 3* RATE" => "PROD NONE", "INJ WATER RATE 86.4" => "INJ WATER NONE",
            "PROD ORAT 43.2 3* RATE" => "PROD NONE", "INJ WATER RATE 43.2" => "INJ WATER NONE"))
        @test isempty(setup_case_from_data_file(no_group).model[:Facility].domain.group_to_wells)
        moved = deepcopy(data)
        moved["SCHEDULE"]["STEPS"][2]["WELSPECS"] = [["P1", "FIELD"]]
        @test_throws "leaves an active allocation group" JutulDarcy.parse_deck_group_steps(moved["SCHEDULE"])
        guide = deepcopy(data)
        guide["SCHEDULE"]["STEPS"][1]["GUIDERAT"] = ["0 OIL 1 0.5 1 1 0 0 YES 0.5"]
        @test_throws "GUIDERAT formulas" JutulDarcy.parse_deck_group_steps(guide["SCHEDULE"])
    end
end
