using Jutul, JutulDarcy, Test, Random, TOML, Statistics
Base.include(JutulDarcy, joinpath(@__DIR__, "cpu_convergence_candidate.jl"))
Random.seed!(1234)
results = []
@testset "Fused CPU convergence reductions" begin
    for phases in (2,3), tagged in (false,true), blackoil in (false,true)
        r = randn(phases, 100_003)
        volume = 1 .+ rand(100_003)
        properties = 1 .+ rand(phases, 100_003)
        if tagged
            volume = map(v -> Jutul.get_ad_entity_scalar(v, phases; tag=Cells()), volume)
            properties = map(v -> Jutul.get_ad_entity_scalar(v, phases; tag=Cells()), properties)
        end
        density = (1000.0,800.0,100.0)
        reference = blackoil ? JutulDarcy.cnv_mb_errors_bo(r,volume,properties,0.25,density,Val(phases)) :
            JutulDarcy.cnv_mb_errors(r,volume,properties,0.25,Val(phases))
        actual = JutulDarcy.cpu_convergence_errors(r,volume,properties,0.25,Val(phases),Val(blackoil),density)
        @test collect(actual[1]) ≈ collect(reference[1])
        @test collect(actual[2]) ≈ collect(reference[2])
        @test JutulDarcy.cpu_convergence_supported(r,volume,properties) == (Threads.nthreads(:default)>1)
        bytes = @allocated JutulDarcy.cpu_convergence_errors(r,volume,properties,0.25,Val(phases),Val(blackoil),density)
        push!(results, Dict("phases"=>phases,"tagged"=>tagged,"blackoil"=>blackoil,"bytes"=>bytes))
    end
    external = fill(Jutul.ForwardDiff.Dual{Nothing}(1.0,1.0),2,20_000)
    @test !JutulDarcy.cpu_convergence_supported(zeros(2,20_000),ones(20_000),external)
    @test_throws DimensionMismatch JutulDarcy.cpu_convergence_errors(zeros(2,3),ones(4),ones(2,3),1.0,Val(2),Val(false))
    nan_residual = fill(NaN,2,20_000)
    nan_result = JutulDarcy.cpu_convergence_errors(nan_residual,ones(20_000),ones(2,20_000),1.0,Val(2),Val(false))
    @test all(isnan,nan_result[1]) && all(isnan,nan_result[2])
end
TOML.print(stdout, Dict("allocation_checks"=>results))
