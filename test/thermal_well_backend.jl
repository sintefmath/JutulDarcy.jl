using Jutul, JutulDarcy, CUDA, Test

@testset "Thermal tracer wells on ka_cuda" begin
    if CUDA.functional()
        CUDA.allowscalar(false)
        path = joinpath(@__DIR__, "data", "thermal", "THERMAL_TRACER_SMALL.DATA")
        case = setup_case_from_data_file(path; dp_max_abs_well = Inf,
            parse_arg = (silent = true,))
        options = (info_level = -1, error_on_incomplete = true,
            nonlinear_tolerance = 1e-7, tol_cnv = 1e-6, tol_mb = 1e-9,
            tol_cnv_well = 1e-7)
        reference = simulate_reservoir(case; options...)
        result = simulate_reservoir(case; mode = :ka_cuda,
            wells_on_device = true, mixed_cross_terms_on_host = true,
            linear_solver_arg = (smoother_type = :ka_dilu,), options...)
        @test length(result.states) == length(case.dt)
        for (actual, expected) in zip(result.states, reference.states)
            @test actual[:Pressure] ≈ expected[:Pressure] rtol = 1e-5
            @test actual[:Temperature] ≈ expected[:Temperature] rtol = 1e-5
            @test all(isfinite, actual[:TracerConcentrations])
            @test actual[:TracerMasses] ≈ expected[:TracerMasses] rtol = 1e-4 atol = 1e-7
        end
    else
        @test_skip CUDA.functional()
    end
end
