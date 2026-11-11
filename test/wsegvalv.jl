using Jutul, JutulDarcy, GeoEnergyIO, Test

struct ConstructorTestSegment{I} end

@testset "Segment valve pressure loss" begin
    JD = JutulDarcy
    valve = SegmentValve(0.9, 2e-4; pipe_length = 0)
    drop(v, q) = segment_pressure_drop(v, 10.0, 0.00065, 0.075, 0.0, q, 1000.0, 0.001)
    @test drop(valve, 2.0) ≈ 2.0^2/(2*1000*(0.9*2e-4)^2)
    @test drop(valve, -2.0) ≈ -drop(valve, 2.0)
    @test drop(valve, 0.0) == 0.0
    @test drop(SegmentValve(0.9, 1e-4; pipe_length = 0), 2.0) ≈ 4*drop(valve, 2.0)
    pipe = SegmentValve(0.9, 2e-4)
    @test drop(pipe, 2.0) ≈ drop(valve, 2.0) + drop(SegmentWellBoreFrictionHB(), 2.0)
    @test drop(SegmentValve(0.9, 2e-4; enabled = false), 2.0) ≈ drop(SegmentWellBoreFrictionHB(), 2.0)
    explicit = SegmentValve(0.9, 2e-4; pipe_length = 12, pipe_diameter = 0.12,
        pipe_roughness = 1e-4, pipe_area = 0.02)
    pipe_loss = segment_pressure_drop(SegmentWellBoreFrictionHB(), 12.0, 1e-4,
        0.06, 0.0, 2.0*(π*0.06^2/0.02), 1000.0, 0.001)
    @test drop(explicit, 2.0) ≈ drop(valve, 2.0) + pipe_loss
    @test isfinite(drop(SegmentValve(0.9, 0; pipe_length = 0), 2.0))
    @test JD.well_segment_is_closed(SegmentValve(0.9, 2e-4; is_open = false))
    @test_throws ArgumentError SegmentValve(0, 2e-4)
    @test_throws ArgumentError SegmentValve(0.9, -1)

    well = MultiSegmentWell([1 2; 2 3], [1, 2], [2, 3];
        segment_models = Any[SegmentWellBoreFrictionHB(), valve])
    @test eltype(well.segment_models) == Union{SegmentWellBoreFrictionHB, SegmentValve}
    for supplied in (Any[valve, valve], Union{SegmentWellBoreFrictionHB, SegmentValve}[valve, valve])
        homogeneous = MultiSegmentWell([1 2; 2 3], [1, 2], [2, 3]; segment_models = supplied)
        @test eltype(homogeneous.segment_models) == SegmentValve
        @test homogeneous.segment_models == supplied
    end
    supplied = Any[ConstructorTestSegment{i}() for i in 1:5]
    many = MultiSegmentWell([1 2 3 4 5; 2 3 4 5 6], [1], [6]; segment_models = supplied)
    @test length(Base.uniontypes(eltype(many.segment_models))) == 5
    @test many.segment_models == supplied
end

@testset "WSEGVALV schedule and simulation" begin
    path = joinpath(@__DIR__, "data", "WSEGVALV.DATA")
    case, data = setup_case_from_data_file(path; include_data = true,
        block_backend = false, parse_arg = (silent = true,))
    well = physical_representation(case.model[:P])
    @test well.neighborship == [2 3; 1 2]
    @test well.perforations.self == [2, 3]
    @test eltype(well.segment_models) == SegmentValve
    @test all(v -> v isa SegmentValve, well.segment_models)
    valve1, valve2 = well.segment_models
    @test valve1.pipe_length == 0
    @test valve2.pipe_length == 10
    @test valve2.pipe_diameter == 0.15
    @test valve2.pipe_roughness == 0.00065
    @test valve2.pipe_area ≈ π*(0.15/2)^2
    @test valve2.maximum_area == valve2.pipe_area
    @test [f[:P].segment_models[1].is_open for f in case.forces] == [true, false, true]
    @test [f[:P].segment_models[1].area for f in case.forces] == [2e-4, 5e-5, 5e-5]
    @test all(f -> f[:P].segment_models[2] == valve2, case.forces)

    result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true,
        timesteps = :none, nonlinear_tolerance = 1e-7)
    @test length(result.states) == 3
    for state in result.result.states
        @test all(isfinite, state[:P][:TotalMassFlux])
        @test all(isfinite, state[:P][:Pressure])
    end
    @test abs(result.result.states[2][:P][:TotalMassFlux][1]) < 1e-10
    @test abs(result.result.states[1][:P][:TotalMassFlux][1]) > 1e-6
    @test abs(result.result.states[3][:P][:TotalMassFlux][1]) > 1e-6
    @test well.segment_models[1].is_open
    @test !haskey(case.state0[:P], :SegmentModels)
    @test all(state -> !haskey(state[:P], :SegmentModels), result.result.states)

    @testset "Merging preserves valve types and forces" begin
        wm = JutulDarcy.WellMerging
        pipe_well = MultiSegmentWell(well.neighborship, well.perforations.reservoir,
            well.perforations.self; name = :Q)
        pipe_domain = DataDomain(pipe_well)
        original_domain = case.model[:P].data_domain
        for (key, (value, entity)) in original_domain.data
            pipe_domain[key, entity] = deepcopy(value)
        end
        merged = wm.merge_well_domains([original_domain, pipe_domain], :SegmentWells)
        merged_well = physical_representation(merged)
        @test eltype(merged_well.segment_models) == Union{SegmentWellBoreFrictionHB, SegmentValve}
        @test count(v -> v isa SegmentValve, merged_well.segment_models) == 2
        forces = Dict(:P => case.forces[2][:P], :Q => (mask = nothing,))
        model = (models = Dict(:SegmentWells => merged),)
        replacements = Dict(:P => :SegmentWells, :Q => :SegmentWells)
        combined = wm.merge_case_forces(forces, model, replacements)[:SegmentWells]
        @test !combined.segment_models[1].is_open
        @test combined.segment_models[2] == valve2
        @test all(v -> v isa SegmentWellBoreFrictionHB, combined.segment_models[3:4])
    end

    @testset "Initially shut valve reopens" begin
        initially_shut = deepcopy(data)
        first(initially_shut["SCHEDULE"]["STEPS"])["WSEGVALV"][1][9] = "SHUT"
        initially_shut["SCHEDULE"]["STEPS"][2]["WSEGVALV"][1][9] = "OPEN"
        shut_case = setup_case_from_data_file(initially_shut; block_backend = false)
        reopened = simulate_reservoir(shut_case; info_level = -1,
            error_on_incomplete = true, timesteps = :none, nonlinear_tolerance = 1e-7)
        @test abs(reopened.result.states[1][:P][:TotalMassFlux][1]) < 1e-10
        @test abs(reopened.result.states[2][:P][:TotalMassFlux][1]) > 1e-6
    end

    @testset "Valves introduced after the first step" begin
        delayed = deepcopy(data)
        steps = delayed["SCHEDULE"]["STEPS"]
        steps[2]["WSEGVALV"] = pop!(steps[1], "WSEGVALV")
        delayed_case = setup_case_from_data_file(delayed; block_backend = false)
        @test eltype(physical_representation(delayed_case.model[:P]).segment_models) == SegmentValve
        @test all(v -> v isa SegmentValve && !v.enabled, delayed_case.forces[1][:P].segment_models)
        @test all(v -> v isa SegmentValve && v.enabled, delayed_case.forces[2][:P].segment_models)
        delayed_result = simulate_reservoir(delayed_case; info_level = -1,
            error_on_incomplete = true, timesteps = :none, nonlinear_tolerance = 1e-7)
        @test length(delayed_result.states) == 3
    end

    @testset "Invalid segment number" begin
        bad = deepcopy(data)
        first(bad["SCHEDULE"]["STEPS"])["WSEGVALV"][1][2] = 99
        @test_throws ArgumentError setup_case_from_data_file(bad)
    end
end
