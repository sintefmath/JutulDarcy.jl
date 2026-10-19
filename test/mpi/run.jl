# Load the packages before starting MPI so that precompilation runs outside
# mpiexec, including on the first run in a fresh environment.
using Test, Jutul, JutulDarcy, MPI, PartitionedArrays, HYPRE

# Use MPI.jl's launcher so that the executable matches the loaded MPI library.
mktempdir() do output_path
    script = joinpath(@__DIR__, "runtests.jl")
    process = run(ignorestatus(`$(MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(Base.active_project()) $script $output_path`))
    success(process) || error("MPI regression test failed (exit code $(process.exitcode)).")
end
