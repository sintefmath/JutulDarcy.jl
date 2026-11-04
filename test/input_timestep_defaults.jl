using Test, Jutul, GeoEnergyIO, JutulDarcy

@testset "Input simulation defaults" begin
    defaults = JutulDarcy.reservoir_simulation_input_defaults
    standard = (
        initial_dt = si_unit(:day),
        max_timestep = si_unit(:year),
        min_timestep = 0.0,
        timestep_max_increase = 10.0,
        timestep_max_decrease = 0.1,
        timestep_decrease_on_cut = 2.0,
        tol_cnv = 1e-3,
        tol_mb = 1e-7,
    )
    @testset "No input settings" begin
        for input in (nothing, missing, (;), Dict(), Dict("SCHEDULE" => Dict()),
                Dict("SCHEDULE" => Dict("STEPS" => [])),
                Dict("SCHEDULE" => Dict("STEPS" => [Dict()])))
            for (key, value) in pairs(standard)
                @test defaults(input, key) == value
            end
        end
    end

    data = parse_data_file(GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA"); silent = true)
    step = first(data["SCHEDULE"]["STEPS"])
    step["TSTEP"] = [86400.0]
    untuned = setup_case_from_data_file(deepcopy(data))
    for (key, value) in pairs(standard)
        @test defaults(untuned, key) == value
    end

    timestep_record = [0.1*86400, 0.25*86400, 0.01*86400, 0.15*86400, 2.0, 0.25, 0.2]
    # Keep supporting cases parsed with the original first-record-only parser.
    step["TUNING"] = timestep_record[1:5]
    legacy = setup_case_from_data_file(deepcopy(data))
    @test defaults(legacy, :initial_dt) == 8640.0
    @test defaults(legacy, :max_timestep) == 21600.0
    @test defaults(legacy, :min_timestep) == 864.0
    @test defaults(legacy, :timestep_max_increase) == 2.0
    @test defaults(legacy, :timestep_max_decrease) == standard.timestep_max_decrease
    @test defaults(legacy, :timestep_decrease_on_cut) == standard.timestep_decrease_on_cut
    @test defaults(legacy, :tol_cnv) == standard.tol_cnv
    @test defaults(legacy, :tol_mb) == standard.tol_mb

    step["TUNING"] = [timestep_record, [0.1, 2e-3, 3e-7]]
    tuned = mktemp() do _, io
        case = redirect_stdout(io) do
            setup_case_from_data_file(data)
        end
        seekstart(io)
        @test !occursin("Keyword TUNING was present, but is not supported", read(io, String))
        return case
    end
    @test defaults(tuned, :tol_cnv) == 2e-3
    @test defaults(tuned.input_data, :tol_mb) == 3e-7

    @testset "Direct simulator setup" begin
        sim, cfg = setup_reservoir_simulator(tuned; info_level = -1)
        selector = first(cfg[:timestep_selectors])
        @test selector.init_abs == 8640.0
        @test selector.max == 21600.0
        @test selector.decrease == 5.0
        @test cfg[:max_timestep] == 21600.0
        @test cfg[:min_timestep] == 864.0
        @test cfg[:timestep_max_increase] == 2.0
        @test cfg[:timestep_max_decrease] == 0.25
        @test Jutul.cut_timestep(sim, cfg, 10000.0, 86400.0, tuned.forces[1], []) == 2000.0
        tol = cfg[:tolerances][:Reservoir][:mass_conservation]
        @test tol.CNV == 2e-3
        @test tol.MB == 3e-7
        for name in keys(get_model_wells(tuned.model))
            well_tol = cfg[:tolerances][name][:mass_conservation]
            @test well_tol.CNV == 10*tol.CNV
        end
        result = simulate_reservoir(tuned; simulator = sim, config = cfg)
        @test first(report_timesteps(result.result.reports; ministeps = true)) == 8640.0
    end

    @testset "Explicit overrides" begin
        for limit in ((max_timestep = 43200.0,), (max_dt = 43200.0,),
                (max_dt = 86400.0, max_timestep = 43200.0),
                (max_dt = 43200.0, max_timestep = 86400.0))
            sim, cfg = setup_reservoir_simulator(tuned; info_level = -1,
                initial_dt = 43200.0, min_timestep = 10.0, timestep_max_increase = 4.0,
                timestep_max_decrease = 0.05, timestep_decrease_on_cut = 4.0,
                tol_cnv = 5e-4, tol_mb = 6e-8, limit...)
            @test cfg[:max_timestep] == 43200.0
            @test cfg[:min_timestep] == 10.0
            @test cfg[:timestep_max_increase] == 4.0
            @test cfg[:timestep_max_decrease] == 0.05
            @test first(cfg[:timestep_selectors]).init_abs == 43200.0
            @test first(cfg[:timestep_selectors]).max == 43200.0
            @test first(cfg[:timestep_selectors]).decrease == 4.0
            @test Jutul.cut_timestep(sim, cfg, 10000.0, 86400.0, tuned.forces[1], []) == 2500.0
            tol = cfg[:tolerances][:Reservoir][:mass_conservation]
            @test tol.CNV == 5e-4
            @test tol.MB == 6e-8
        end
    end

    @testset "Simulation setup" begin
        result = simulate_reservoir(tuned; info_level = -1)
        cfg = result.extra[:config]
        @test cfg[:max_timestep] == 21600.0
        @test cfg[:min_timestep] == 864.0
        @test cfg[:timestep_max_increase] == 2.0
        @test cfg[:timestep_max_decrease] == 0.25
        @test first(cfg[:timestep_selectors]).decrease == 5.0
        @test cfg[:tolerances][:Reservoir][:mass_conservation].CNV == 2e-3
        @test cfg[:tolerances][:Reservoir][:mass_conservation].MB == 3e-7
        @test first(report_timesteps(result.result.reports; ministeps = true)) == 8640.0
        result = simulate_reservoir(tuned; info_level = -1, max_dt = 43200.0,
            initial_dt = 43200.0, timestep_max_increase = 4.0)
        @test result.extra[:config][:max_timestep] == 43200.0
        @test result.extra[:config][:timestep_max_increase] == 4.0
        @test first(report_timesteps(result.result.reports; ministeps = true)) == 43200.0
    end

    @testset "Untuned setup" begin
        _, cfg = setup_reservoir_simulator(untuned; info_level = -1)
        @test cfg[:max_timestep] == standard.max_timestep
        @test cfg[:min_timestep] == standard.min_timestep
        @test cfg[:timestep_max_increase] == standard.timestep_max_increase
        @test cfg[:timestep_max_decrease] == standard.timestep_max_decrease
        @test first(cfg[:timestep_selectors]).init_abs == standard.initial_dt
        @test first(cfg[:timestep_selectors]).decrease == standard.timestep_decrease_on_cut
        @test cfg[:tolerances][:Reservoir][:mass_conservation].CNV == standard.tol_cnv
        @test cfg[:tolerances][:Reservoir][:mass_conservation].MB == standard.tol_mb
    end

    @testset "Later TUNING changes" begin
        later = deepcopy(data)
        push!(later["SCHEDULE"]["STEPS"], Dict("TUNING" => [2 .* timestep_record, [0.1, 1e-2, 1e-6]]))
        for key in keys(standard)
            @test defaults(later, key) == defaults(tuned, key)
        end
        later_case = JutulCase(tuned.model, tuned.dt, tuned.forces;
            state0 = tuned.state0, parameters = tuned.parameters, input_data = later)
        mktemp() do _, io
            redirect_stdout(io) do
                setup_reservoir_simulator(later_case; info_level = -1)
            end
            seekstart(io)
            @test count("later TUNING changes are ignored.", read(io, String)) == 1
        end
        pop!(first(later["SCHEDULE"]["STEPS"]), "TUNING")
        for (key, value) in pairs(standard)
            @test defaults(later, key) == value
        end
    end
end
