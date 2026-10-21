# Usage: julia compare_states.jl reference.state candidate.state [candidate.state...]
using Serialization, LinearAlgebra, TOML

reference = deserialize(first(ARGS))
comparisons = Any[]
for path in ARGS[2:end]
    candidate = deserialize(path)
    @assert size(reference.pressure) == size(candidate.pressure)
    @assert size(reference.saturations) == size(candidate.saturations)
    @assert all(isfinite, candidate.pressure) && all(isfinite, candidate.saturations)
    dp = candidate.pressure .- reference.pressure
    ds = candidate.saturations .- reference.saturations
    row = Dict("reference" => basename(first(ARGS)), "candidate" => basename(path),
        "pressure_max_absolute" => maximum(abs, dp),
        "pressure_relative_norm" => norm(dp)/norm(reference.pressure),
        "saturation_max_absolute" => maximum(abs, ds),
        "saturation_relative_norm" => norm(ds)/norm(reference.saturations),
        "bitwise_equal" => isequal(reference, candidate))
    push!(comparisons, row)
    println(row)
end
if haskey(ENV, "PERF_COMPARISON_OUTPUT")
    open(ENV["PERF_COMPARISON_OUTPUT"], "w") do io
        TOML.print(io, Dict("comparisons" => comparisons))
    end
end
