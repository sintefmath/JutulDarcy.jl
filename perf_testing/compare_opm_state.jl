# Usage: julia compare_opm_state.jl reference.state phase_metadata.toml reference.cells CASE.UNRST
# Align OPM fields to the Jutul mesh before using the existing whole-state check.
using Serialization, TOML
reference_path, phase_metadata_path, cell_map_path, restart_path = ARGS
prefix = joinpath(dirname(restart_path), "final_aligned")
run(`/usr/bin/python3 $(joinpath(@__DIR__, "extract_opm_fields.py")) $restart_path $prefix`)
metadata = TOML.parsefile(prefix*".fields.toml")
reference_metadata = TOML.parsefile(phase_metadata_path)
full_metadata = TOML.parsefile(splitext(reference_path)[1]*".toml")
@assert metadata["report"] == full_metadata["report_steps"] "Restart does not contain the requested final report"
@assert isapprox(metadata["elapsed_days"],sum(full_metadata["report_days"]); rtol=1e-13,atol=1e-7) "Restart physical time differs from the reference"
reference_cells = deserialize(cell_map_path)
n = metadata["cells"]
opm_cells = read!(prefix*".cells.bin", Vector{Int64}(undef,n))
@assert sort(reference_cells) == sort(opm_cells) "Active-cell identities differ between the simulators"
positions = zeros(Int, maximum(opm_cells))
positions[opm_cells] .= 1:n
permutation = positions[reference_cells]
@assert all(>(0), permutation)
phase_indices = [only(findall(==(name), metadata["phase_names"])) for name in reference_metadata["phase_names"]]
pressure = read!(prefix*".pressure.bin", Vector{Float64}(undef,n))[permutation]
saturations = reshape(read!(prefix*".saturations.bin", Vector{Float64}(undef,3n)),3,n)[phase_indices,permutation]
serialize(prefix*".state", (; pressure,saturations))
ENV["PERF_COMPARISON_OUTPUT"] = prefix*".comparison.toml"
empty!(ARGS)
append!(ARGS,[reference_path,prefix*".state"])
include(joinpath(@__DIR__,"compare_states.jl"))
