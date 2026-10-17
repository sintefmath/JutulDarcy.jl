using Test, Jutul, GeoEnergyIO, JutulDarcy

@testset "Input timestep defaults" begin
    defaults = JutulDarcy.reservoir_simulation_input_defaults
    @test defaults(nothing) == (;)
    @test defaults(missing) == (;)
    @test defaults(Dict()) == (;)
    data = parse_data_file(GeoEnergyIO.test_input_file_path("SPE1","SPE1.DATA");silent=true)
    first(data["SCHEDULE"]["STEPS"])["TSTEP"] = [86400.0]
    untuned = setup_case_from_data_file(deepcopy(data))
    first(data["SCHEDULE"]["STEPS"])["TUNING"] = [0.1*86400,0.25*86400,0.01*86400,0.15*86400,2.0]
    tuned = setup_case_from_data_file(data)
    @test defaults(tuned.input_data) == (initial_dt=8640.,max_timestep=21600.,min_timestep=864.,timestep_max_increase=2.)
    result = simulate_reservoir(tuned;info_level=-1)
    cfg = result.extra[:config]
    @test cfg[:max_timestep] == 21600.
    @test cfg[:min_timestep] == 864.
    @test cfg[:timestep_max_increase] == 2.
    @test first(report_timesteps(result.result.reports;ministeps=true)) == 8640.
    result = simulate_reservoir(tuned;info_level=-1,max_dt=43200.,initial_dt=43200.,timestep_max_increase=4.)
    @test result.extra[:config][:max_timestep] == 43200.
    @test result.extra[:config][:timestep_max_increase] == 4.
    @test first(report_timesteps(result.result.reports;ministeps=true)) == 43200.
    result = simulate_reservoir(untuned;info_level=-1)
    @test result.extra[:config][:max_timestep] == si_unit(:year)
    @test result.extra[:config][:min_timestep] == 0.
    @test result.extra[:config][:timestep_max_increase] == 10.
end
