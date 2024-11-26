using Test, JutulDarcy
@testset "Discretizations" begin
    pth = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
    subs = 1:2
    for kgrad in [:tpfa, :tpfa_test, :avgmpfa, :ntpfa]
        for upw in [:spu, :weno]
            @testset "$kgrad-$upw" begin
                case = setup_case_from_data_file(pth,
                    kgrad = kgrad,
                    upwind = upw,
                )
                _, states = simulate_reservoir(case[subs],
                    info_level = -1,
                    failure_cuts_timestep = false
                )
                @test length(states) == length(subs)
            end
        end
    end
end
