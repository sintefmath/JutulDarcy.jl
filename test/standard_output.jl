using Dates, GeoEnergyIO, Jutul, JutulDarcy, Test

@testset "Standard simulation output" begin
    deck_path = joinpath(@__DIR__, "data", "thermal", "HEAT_SMALL.DATA")
    data = parse_data_file(deck_path; silent = true)
    dims = Tuple(data["GRID"]["cartDims"])
    actnum = trues(dims)
    actnum[2, 1, 1] = false
    data["GRID"]["ACTNUM"] = actnum
    case = setup_case_from_data_file(data; verbose = false)[1:2]
    result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true)
    nc = count(actnum)

    mktempdir() do dir
        skipped = write_standard_output(result, case, joinpath(dir, "SKIPPED");
            summary = false, restart = false, egrid = false)
        @test skipped == (summary = nothing, restart = nothing, egrid = nothing)
        @test isempty(readdir(dir))

        path = joinpath(dir, "DECK")
        paths = write_standard_output(result, case, path)
        grid = read_egrid(paths.egrid; units = :native)
        @test Tuple(grid["cartDims"]) == dims
        @test grid["ACTNUM"] == actnum
        @test first(grid["GRIDUNIT"]) == "FEET"
        @test maximum(grid["COORD"][1:3:end]) ≈ 100.0

        summary = read_summary(path; format = :structured)
        @test summary["UNIT_SYSTEM"] == "FIELD"
        @test summary["DIMENS"] == dims
        @test summary["TIME"].start_date == DateTime(1994, 1, 1)
        @test summary["TIME"].seconds ≈ result.time rtol = 1e-6
        @test summary["VALUES"]["WELLS"]["P1"]["WBHP"] ≈ result.wells[:P1][:bhp]/si_unit(:psi) rtol = 1e-6

        steps = read_restart(paths.restart; units = :si)
        @test length(steps) == length(result.states)
        for (i, step) in enumerate(steps)
            state = result.states[i]
            @test step["report_step"] == i
            @test step["INTEHEAD"][3] == 2
            @test step["INTEHEAD"][9:12] == [dims..., nc]
            @test step["days"] ≈ result.time[i]/si_unit(:day)
            @test step["PRESSURE"] ≈ state[:Pressure] rtol = 1e-6
            @test step["SWAT"] ≈ state[:Saturations][1, :] rtol = 1e-6
            @test step["RS"] ≈ state[:Rs] rtol = 1e-6
            @test step["TEMP"] ≈ state[:Temperature] rtol = 1e-6
        end

        lab_path = joinpath(dir, "LAB")
        lab = write_standard_output(result, case, lab_path; units = :lab, unified = false, restart = false)
        @test isnothing(lab.restart)
        @test !isfile(lab_path * ".UNRST")
        @test length(lab.summary.data) == length(result.time)
        @test first(read_egrid(lab.egrid; units = :native)["GRIDUNIT"]) == "CM"
        lab_summary = read_summary(lab_path; format = :structured)
        @test lab_summary["UNIT_SYSTEM"] == "LAB"
        @test lab_summary["VALUES"]["WELLS"]["I1"]["WWIR"] ≈ result.wells[:I1][:wrat]*si_unit(:hour)/1e-6 rtol = 1e-6

        synthetic_case = JutulCase(case.model, case.dt, case.forces;
            state0 = case.state0, parameters = case.parameters, start_date = DateTime(2020, 5, 1))
        synthetic_path = joinpath(dir, "SYNTHETIC")
        synthetic = write_standard_output(result, synthetic_case, synthetic_path; unified = false)
        @test read_egrid(synthetic.egrid)["cartDims"] == [nc, 1, 1]
        @test length(synthetic.restart) == length(steps)
        @test read_restart(synthetic_path)[1]["INTEHEAD"][3] == 1
        @test read_restart(synthetic_path)[1]["INTEHEAD"][9:12] == [nc, 1, 1, nc]
        synthetic_summary = read_summary(synthetic_path; format = :structured)
        @test synthetic_summary["DIMENS"] == (nc, 1, 1)
        @test synthetic_summary["TIME"].start_date == synthetic_case.start_date
    end
end

@testset "Black-oil restart ratios and phase properties" begin
    jd = JutulDarcy
    pvtg = jd.PVTG(jd.PVTGTable([[20.0, 21.0], [100.0, 105.0]],
        [[0.0, 0.001], [0.0, 0.01]], [[1e-5, 1e-5], [2e-5, 2e-5]],
        [1e6, 2e7], [0.001, 0.01]))
    domain = reservoir_domain(CartesianMesh((3, 1, 1)))
    model, parameters = jd.setup_reservoir_model_from_blackoil_tables(domain;
        pvto = jd.spe1_pvto(), pvtg, pvtw = jd.spe1_pvtw(),
        reference_densities = (1000.0, 800.0, 1.0), extra_out = true, extra_outputs = true)
    saturations = repeat([0.2, 0.5, 0.3], 1, 3)
    state0 = setup_reservoir_state(model; Pressure = fill(1e7, 3), Saturations = saturations,
        Rs = zeros(3), Rv = zeros(3))
    case = JutulCase(model, [1.0]; state0, parameters)
    result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true)
    state = only(result.states)
    @test all(>(0), state[:Rs])
    @test all(>(0), state[:Rv])
    gas_per_oil = 1000*si_unit(:feet)^3/si_unit(:stb)
    density = si_unit(:pound)/si_unit(:feet)^3

    mktempdir() do dir
        for (units, rs_factor, pressure, rho_factor) in ((:field, gas_per_oil, si_unit(:psi), density),
                (:metric, 1.0, si_unit(:bar), 1.0), (:lab, 1.0, si_unit(:atm), 1000.0))
            path = write_standard_restart(result, case, joinpath(dir, string(units)); units)
            step = only(read_restart(path))
            @test step["PRESSURE"] ≈ state[:Pressure]/pressure rtol = 1e-6
            @test step["RS"] ≈ state[:Rs]/rs_factor rtol = 1e-6
            @test step["RV"] ≈ state[:Rv]*rs_factor rtol = 1e-6
            @test step["OIL_DEN"] ≈ state[:PhaseMassDensities][2, :]/rho_factor rtol = 1e-6
            @test step["GAS_VISC"] ≈ state[:PhaseViscosities][3, :]/1e-3 rtol = 1e-6
            @test step["1OVERBG"] ≈ state[:ShrinkageFactors][3, :]/rs_factor rtol = 1e-6
            @test step["WATKR"] ≈ state[:RelativePermeabilities][1, :] rtol = 1e-6
        end

        reduced = result.result.states[1][:Reservoir]
        expected = Dict(k => pop!(reduced, k) for k in (:Saturations, :Rs, :Rv))
        paths = write_standard_output(result, case, joinpath(dir, "REDUCED"); units = :metric)
        step = only(read_restart(paths.restart))
        @test step["SOIL"] ≈ expected[:Saturations][2, :] rtol = 1e-6
        @test step["RS"] ≈ expected[:Rs] rtol = 1e-6
        @test step["RV"] ≈ expected[:Rv] rtol = 1e-6
        summary = read_summary(paths.summary.spec; format = :structured)
        @test summary["VALUES"]["FIELD"]["FOIPR"] ≈ result.summary["VALUES"]["FIELD"]["FOIPR"] rtol = 1e-6
        @test !haskey(reduced, :Rs)
    end
end

@testset "Compositional restart mole fractions and isothermal temperature" begin
    deck = GeoEnergyIO.test_input_file_path("SIMPLE_COMP", "SIMPLE_COMP.DATA")
    # Use a short row with the same EOS, rock and fluid properties.
    data = mktempdir() do dir
        path = joinpath(dir, "SMALL_COMP.DATA")
        input = replace(read(deck, String), "1000*" => "8*", "1000 1 1" => "8 1 1", "PROD * 1000 1 0" => "PROD * 8 1 0")
        write(path, input)
        parse_data_file(path; silent = true)
    end
    dims = Tuple(data["GRID"]["cartDims"])
    @test dims == (8, 1, 1)
    case = setup_case_from_data_file(data; extra_outputs = true)[1:2]
    result = simulate_reservoir(case; info_level = -1, error_on_incomplete = true)
    @test case.parameters[:Reservoir][:Temperature] ≈ fill(423.15, prod(dims))
    eos = reservoir_model(case).system.equation_of_state
    mw = JutulDarcy.MultiComponentFlash.molar_masses(eos)

    mktempdir() do dir
        for (units, temperature) in ((:field, 302.0), (:metric, 150.0), (:lab, 150.0))
            path = write_standard_restart(result, case, joinpath(dir, string(units)); units)
            for (step, state) in zip(read_restart(path), result.states)
                @test all(x -> isapprox(x, temperature; rtol = 1e-6), step["TEMP"])
                @test !haskey(step, "RS") && !haskey(step, "RV")
                for (name, key) in (("ZMF", :OverallMoleFractions), ("XMF", :LiquidMassFractions), ("YMF", :VaporMassFractions))
                    expected = state[key]
                    if name != "ZMF"
                        expected = expected ./ collect(mw)
                        expected = expected ./ sum(expected; dims = 1)
                    end
                    for i in eachindex(mw)
                        @test step["$name$i"] ≈ expected[i, :] rtol = 1e-6
                    end
                    @test sum(step["$name$i"] for i in eachindex(mw)) ≈ ones(prod(dims)) rtol = 1e-6
                end
            end
        end
    end
end
