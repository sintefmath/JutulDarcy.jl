# MPI support

There are a few adjustments needed before a script can be run in MPI.

## Setting up the environment

You will have to set up an environment with the following packages under Julia 1.9+:
`PartitionedArrays`, `MPI`, and `JutulDarcy`. `HYPRE` is automatically installed. This is generally the best performing solver setup available, even if you are working in a shared memory environment.

### Writing the script

Write your script as usual. The following options must then be set:

- [`setup_reservoir_model`](@ref) should have the extra keyword argument `split_wells=true`. We also recommend keeping the default `backend=:csr` for the best performance.
- [`simulate_reservoir`](@ref) or [`setup_reservoir_simulator`](@ref) should get the optional argument `mode = :mpi`

You must then run the file using the appropriate `mpiexec` as described in the MPI.jl documentation. Specialized functions will be called by `simulate_reservoir` when this is the case. We document them here, even if we recommend using the high level version of this interface:

```@docs
JutulDarcy.simulate_reservoir_parray
```

## Checklist for running in MPI

- Install and load the following packages at the top of your script: `PartitionedArrays, MPI, HYPRE`
- (Recommended): Put `MPI.Init()` at the top of your script
- Make sure that `split_wells = true` is set in your model setup
- Set `output_path` when running the simulation (otherwise the results will not be stored anywhere)
- Set `mode=:mpi` when running the simulation.

A few useful functions:

- `MPI.install_mpiexecjl()` installs MPI that works "out of the box" with your current Julia setup
- `MPI.mpiexec()` gives you the path to the executable and all environment variables

A typical command to launch a MPI script from within Julia:

```julia
n = 5 # = 5 processes
script_to_run = "my_script.jl"
run(`$(mpiexec()) -n $n $(Base.julia_cmd()) --project=$(Base.active_project()) $script_to_run`)
```

Adding threads to the command will make JutulDarcy use both threads and processes

:warning: **Running a script in MPI means that all parts of the script will run on each process!** If you want to do data analysis you will have to either wrap your code in `if MPI.Comm_rank(MPI.COMM_WORLD) == 0` or do the data analysis in serial (recommended).

## Limitations of running in MPI

MPI can be cumbersome to use when compared to a standard Julia script, and the
current implementation relies on the model being set up on each processor before
subdivision. This can be quite memory intensive during startup.

You should be familiar with the MPI programming model to use this feature. See
[MPI.jl](https://juliaparallel.org/MPI.jl/stable/) for more details, and how MPI
is handled in Julia specifically.

For larger models, compiling the [Standalone reservoir simulator](@ref) is
highly recommended.

!!! note
    MPI consolidates results by writing files to disk. Unless you have a plan to work with the distributed states in-memory returned by the `simulate!` call, it is best to specify a `output_path` optional argument to [`setup_reservoir_simulator`](@ref). After the simulation, that folder will contain output just as if you had run the case in serial.
