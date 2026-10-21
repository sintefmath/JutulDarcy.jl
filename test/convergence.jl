using Jutul, JutulDarcy, Test, Random, TOML, Statistics

function serial_convergence_reference(r, volume, properties, dt, density, ::Val{N}, blackoil) where N
    total = sum(Jutul.value, volume)
    if blackoil
        average = ntuple(p -> sum(v -> inv(Jutul.value(v)), view(properties,p,:))/length(volume), Val(N))
        cnv = ntuple(p -> average[p]*dt*mapreduce((v,pv) -> abs(v)/Jutul.value(pv), max, view(r,p,:), volume)/density[p], Val(N))
        mb = ntuple(p -> average[p]*dt*abs(sum(view(r,p,:)))/(density[p]*total), Val(N))
    else
        cnv = ntuple(p -> JutulDarcy.cnv_error_for_phase(r,volume,properties,dt,p), Val(N))
        mb = ntuple(p -> JutulDarcy.mb_error_for_phase(r,total,properties,dt,p), Val(N))
    end
    return cnv, mb
end

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
        reference = serial_convergence_reference(r,volume,properties,0.25,density,Val(phases),blackoil)
        actual = JutulDarcy.cpu_convergence_errors(r,volume,properties,0.25,Val(phases),Val(blackoil),density)
        @test collect(actual[1]) ≈ collect(reference[1])
        @test collect(actual[2]) ≈ collect(reference[2])
        @test JutulDarcy.cpu_convergence_supported(r,volume,properties) == (Threads.nthreads(:default)>1)
        production = () -> blackoil ? JutulDarcy.cnv_mb_errors_bo(r,volume,properties,0.25,density,Val(phases)) :
            JutulDarcy.cnv_mb_errors(r,volume,properties,0.25,Val(phases))
        @test collect(production()[1]) ≈ collect(reference[1])
        @test collect(production()[2]) ≈ collect(reference[2])
        if Threads.nthreads(:default)>1
            @test (@allocated production()) <= 4096
        end
        bytes = @allocated JutulDarcy.cpu_convergence_errors(r,volume,properties,0.25,Val(phases),Val(blackoil),density)
        push!(results, Dict("phases"=>phases,"tagged"=>tagged,"blackoil"=>blackoil,"bytes"=>bytes))
    end
    external = fill(Jutul.ForwardDiff.Dual{Nothing}(1.0,1.0),2,20_000)
    @test !JutulDarcy.cpu_convergence_supported(zeros(2,20_000),ones(20_000),external)
    @test !JutulDarcy.cpu_convergence_supported(zeros(Float32,2,20_000),ones(20_000),ones(2,20_000))
    @test !JutulDarcy.cpu_convergence_supported(zeros(2,20_000),ones(20_000),ones(2,20_000),Jutul.ForwardDiff.Dual{Nothing}(1.0,1.0))
    @test_throws DimensionMismatch JutulDarcy.cpu_convergence_errors(zeros(2,3),ones(4),ones(2,3),1.0,Val(2),Val(false))
    nan_residual = fill(NaN,2,20_000)
    nan_result = JutulDarcy.cpu_convergence_errors(nan_residual,ones(20_000),ones(2,20_000),1.0,Val(2),Val(false))
    @test all(isnan,nan_result[1]) && all(isnan,nan_result[2])
    negative_volume = fill(-1.0,20_000)
    negative_reference = serial_convergence_reference(ones(2,20_000),negative_volume,ones(2,20_000),1.0,nothing,Val(2),false)
    negative_actual = JutulDarcy.cpu_convergence_errors(ones(2,20_000),negative_volume,ones(2,20_000),1.0,Val(2),Val(false))
    @test negative_reference == negative_actual
end
TOML.print(stdout, Dict("allocation_checks"=>results))
