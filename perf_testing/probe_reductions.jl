# Focused host convergence probe. Numerical reference is the application code.
using Jutul, JutulDarcy, Polyester, Statistics, TOML, Random
using Jutul: ForwardDiff
Random.seed!(1234)

function fused_phase(r, pore_volume, properties, dt, phase, blackoil)
    maximum_residual = 0.0
    residual_sum = 0.0
    property_sum = 0.0
    @batch minbatch=1000 reduction=((max, maximum_residual), (+, residual_sum), (+, property_sum)) for cell in eachindex(pore_volume)
        @inbounds begin
            residual = Jutul.value(r[phase, cell])
            volume = Jutul.value(pore_volume[cell])
            property = Jutul.value(properties[phase, cell])
        end
        residual_sum += residual
        if blackoil
            property_sum += inv(property)
            maximum_residual = max(maximum_residual, abs(residual)/volume)
        else
            property_sum += abs(property)
            maximum_residual = max(maximum_residual, dt*abs(residual)/(property*volume))
        end
    end
    return maximum_residual, residual_sum, property_sum/length(pore_volume)
end

function batched_volume(pore_volume)
    total_volume = 0.0
    @batch minbatch=1000 reduction=((+, total_volume),) for cell in eachindex(pore_volume)
        @inbounds total_volume += Jutul.value(pore_volume[cell])
    end
    return total_volume
end

function fused_errors(r, pore_volume, properties, dt, ::Val{N}, blackoil) where N
    total_volume = batched_volume(pore_volume)
    metrics = ntuple(phase -> fused_phase(r, pore_volume, properties, dt, phase, blackoil), Val(N))
    reference_density = (1000.0, 800.0, 100.0)
    cnv = ntuple(Val(N)) do phase
        maximum_residual, residual_sum, property_mean = metrics[phase]
        blackoil ? property_mean*dt*maximum_residual/reference_density[phase] : maximum_residual
    end
    mb = ntuple(Val(N)) do phase
        maximum_residual, residual_sum, property_mean = metrics[phase]
        blackoil ? property_mean*dt*abs(residual_sum)/(reference_density[phase]*total_volume) :
            (dt/total_volume)*abs(residual_sum)/property_mean
    end
    return cnv, mb
end

function measurement(f; repetitions=50)
    f()
    samples = [@timed(for _ in 1:repetitions; f(); end) for _ in 1:5]
    return Dict("seconds" => median(s.time for s in samples)/repetitions,
        "bytes" => median(s.bytes for s in samples)/repetitions)
end

results = []
for n in (177724, 385500), phases in (2, 3), ad in (false, true)
    r = randn(phases, n)
    volume = 1 .+ rand(n)
    properties = 1 .+ rand(phases, n)
    if ad
        properties = map(v -> Jutul.get_ad_entity_scalar(v, phases; tag=Cells()), properties)
    end
    blackoil = phases == 3
    reference = blackoil ? () -> JutulDarcy.cnv_mb_errors_bo(r, volume, properties, 1.0, (1000.0, 800.0, 100.0), Val(phases)) :
        () -> JutulDarcy.cnv_mb_errors(r, volume, properties, 1.0, Val(phases))
    fused = () -> fused_errors(r, volume, properties, 1.0, Val(phases), blackoil)
    a, b = reference(), fused()
    @assert all(isapprox.(a[1], b[1]; rtol=1e-12))
    @assert all(isapprox.(a[2], b[2]; rtol=1e-10, atol=1e-15)) (a, b)
    row = Dict("cells" => n, "phases" => phases, "ad" => ad,
        "threads" => Threads.nthreads(), "original" => measurement(reference), "fused" => measurement(fused))
    push!(results, row)
    println(row)
    flush(stdout)
end
out = joinpath(@__DIR__, "results", "reduction_t$(Threads.nthreads()).toml")
open(out, "w") do io
    TOML.print(io, Dict("results" => results))
end
