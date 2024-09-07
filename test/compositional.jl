using JutulDarcy
using Jutul
using Test

@testset "SIMPLE_COMP validation" begin
    dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("SIMPLE_COMP")
    data_path = joinpath(dpth, "SIMPLE_COMP.DATA")
    ws, states = result;
    ref_path = joinpath(dpth, "reference.jld2")
    ref = Jutul.JLD2.load(ref_path)
    states_cmp = ref["e300"]
    for fast_flash in (false, true)
        case = setup_case_from_data_file(data_path, fast_flash = fast_flash)
        result = simulate_reservoir(case, info_level = -1)
        @testset "OverallMoleFractions" begin
            for i in 1:500
                rstate = states_cmp[i]
                state = states[i]
                z0 = rstate[:OverallMoleFractions]
                z = state[:OverallMoleFractions]
                @test norm(z-z0)/norm(z) < 0.05
            end
        end
        @testset "Pressure" begin
            for i in 1:500
                rstate = states_cmp[i]
                state = states[i]
                p0 = rstate[:Pressure]
                p = state[:Pressure]
                @test norm(p-p0)/norm(p) < 0.01
            end
        end

        @testset "Saturations" begin
            for i in 1:500
                rstate = states_cmp[i]
                state = states[i]
                s0 = rstate[:Saturations][2, :]
                s = state[:Saturations][2, :]
                @test norm(s-s0)/norm(s) < 0.05
            end
        end
    end
end
