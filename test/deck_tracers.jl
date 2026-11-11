using Jutul, JutulDarcy, GeoEnergyIO, Test

@testset "Drogon passive tracer unit labels" begin
    path = joinpath(@__DIR__, "data", "thermal", "DROGON_TRACERS_SMALL.DATA")
    data = parse_data_file(path; silent = true)
    @test [rec[3] for rec in data["PROPS"]["TRACER"]] == ["g", "g"]
    case = setup_case_from_data_file(deepcopy(data))
    rmodel = reservoir_model(case.model)
    flux = rmodel.equations[:tracers].flux_type
    @test flux.names == ("WT1", "WT2")
    @test all(t -> t isa SinglePhaseTracer, flux.tracers)
    @test number_of_tracers(case.model) == 2
    @test all(iszero, case.state0[:Reservoir][:TracerConcentrations])
    @test [f[:Facility].control[:I].tracers for f in case.forces] ==
        [[3.125, 0.0], [3.125, 0.0], [0.0, 3.125], [0.0, 0.0]]

    # A unit label does not rescale the passive tracer's concentrations.
    default_units = deepcopy(data)
    for rec in default_units["PROPS"]["TRACER"]
        rec[3] = ""
    end
    reference = setup_case_from_data_file(default_units)
    options = (info_level = -1, error_on_incomplete = true, timesteps = :none,
        nonlinear_tolerance = 1e-7, tol_cnv = 1e-6, tol_mb = 1e-9)
    result = simulate_reservoir(case; options...)
    expected = simulate_reservoir(reference; options...)
    @test length(result.states) == length(expected.states) == length(case.dt)
    for (state, ref) in zip(result.states, expected.states)
        @test state[:Pressure] ≈ ref[:Pressure]
        @test state[:TracerConcentrations] ≈ ref[:TracerConcentrations]
        @test all(c -> isfinite(c) && -1e-8 <= c <= 3.125 + 1e-7,
            state[:TracerConcentrations])
    end
    @test all(iszero, result.states[2][:TracerConcentrations][2, :])
    @test all(i -> maximum(result.states[end][:TracerConcentrations][i, :]) > 0, 1:2)

    @testset "Nonzero initial concentrations" begin
        initial = deepcopy(data)
        initial["SOLUTION"]["TVDPFWT1"] = [1000.0 0.0; 1010.0 2.0]
        initial["SOLUTION"]["TBLKFWT2"] = fill(0.75, 5, 5, 1)
        initialized = setup_case_from_data_file(initial)
        concentrations = initialized.state0[:Reservoir][:TracerConcentrations]
        @test concentrations[1, :] ≈ ones(25)
        @test concentrations[2, :] == fill(0.75, 25)
    end

    @testset "Unsupported tracer physics" begin
        for (index, value) in ((4, "OIL"), (5, 1), (6, "ADS"))
            unsupported = deepcopy(data)
            unsupported["PROPS"]["TRACER"][1][index] = value
            @test_throws AssertionError setup_case_from_data_file(unsupported)
        end
    end
end
