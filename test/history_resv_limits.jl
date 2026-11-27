using Test, Jutul, GeoEnergyIO, JutulDarcy

@testset "History RESV limits persist across ministeps" begin
    data = parse_data_file(GeoEnergyIO.test_input_file_path("SPE1","SPE1.DATA");silent=true)
    step = first(data["SCHEDULE"]["STEPS"])
    step["TSTEP"] = [86400.]
    prod = only(pop!(step,"WCONPROD"))
    step["WCONHIST"] = [[prod[1],prod[2],"RESV",prod[4],0.,0.,0,0,0.,prod[9]]]
    case = setup_case_from_data_file(data)
    # SPE1 has DISGAS and no VAPOIL, so the state has no Rv field.
    avg = JutulDarcy.setup_average_resv_state(reservoir_model(case),
        (FluidVolume=[1.],Pressure=[2e7],ImmiscibleSaturation=[0.2],Rs=[10.]))
    @test avg.rv == 0.
    seen = Bool[]
    function check_limits(converged,report,storage,model,dt,forces,cfg,it)
        limits = storage[:Facility].state.WellGroupConfiguration.limits[:PROD]
        push!(seen,!isnothing(limits) && haskey(limits,:resv))
        return converged
    end
    result = simulate_reservoir(case;info_level=-1,initial_dt=21600.,max_timestep=21600.,post_iteration_hook=check_limits)
    @test length(report_timesteps(result.result.reports;ministeps=true)) >= 4
    @test !isempty(seen)
    @test all(seen)
end

@testset "History RESV without Rs" begin
    data = parse_data_file(GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA"); silent = true)
    pop!(data["RUNSPEC"], "DISGAS")
    data["RUNSPEC"]["VAPOIL"] = true
    props = data["PROPS"]
    pop!(props, "PVTO")
    pop!(props, "PVDG")
    props["PVDO"] = [[1e5 1.0 1e-3; 1e8 1.0 1e-3]]
    props["PVTG"] = [Dict("pos" => [1, 3, 5], "key" => [1e5, 1e8],
        "data" => [0.0 0.01 1e-5; 0.001 0.01 1e-5; 0.0 0.002 1e-5; 0.002 0.002 1e-5])]
    step = first(data["SCHEDULE"]["STEPS"])
    step["TSTEP"] = [86400.]
    prod = only(pop!(step, "WCONPROD"))
    step["WCONHIST"] = [[prod[1], prod[2], "RESV", prod[4], 0., 0., 0, 0, 0., prod[9]]]
    case = setup_case_from_data_file(data)
    @test !haskey(case.state0[:Reservoir], :Rs)
    avg = JutulDarcy.setup_average_resv_state(reservoir_model(case),
        (FluidVolume = [1., 3.], Pressure = [1e7, 3e7],
        ImmiscibleSaturation = [0.2, 0.2], Rv = [0.0001, 0.0003]))
    @test avg.rs == 0.
    @test avg.rv ≈ 0.00025
    @test avg.p ≈ 2.5e7
    result = simulate_reservoir(case; info_level = -1, initial_dt = 21600., max_timestep = 21600., error_on_incomplete = true)
    @test length(report_timesteps(result.result.reports; ministeps = true)) >= 4
    @test length(result.states) == 1
    @test !haskey(only(result.states), :Rs)
    @test all(isfinite, only(result.states)[:Pressure])
end
