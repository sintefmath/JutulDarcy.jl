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
