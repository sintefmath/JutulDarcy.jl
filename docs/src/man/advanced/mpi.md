# Multi-threading and MPI support

JutulDarcy can use threads by default, but advanced options can improve performance significantly for larger models.

## Overview of parallel support

There are two main ways of exploiting multiple cores in Jutul/JutulDarcy: Threads are automatically used for assembly and can be used for parts of the linear solve. If you require the best performance, you have to go to MPI where the linear solvers can use a parallel [BoomerAMG preconditioner](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) via [HYPRE.jl](https://github.com/fredrikekre/HYPRE.jl). In addition, there is experimental GPU support described in [GPU support](@ref).

### MPI parallelization

MPI parallelizes all aspects of the solver using domain decomposition and allows a simulation to be divided between multiple nodes in e.g. a supercomputer. It is significantly more cumbersome to use than standard simulations as the program must be launched in MPI mode. This is typically a non-interactive process where you launch your MPI processes and once they complete the simulation the result is available on disk. The MPI parallel option uses a combination of MPI.jl, PartitionedArrays.jl and HYPRE.jl.

### Thread parallelization

JutulDarcy also supports threads. By default, this only parallelizes property evaluations and assembly of the linear system. For many problems, the linear solve is the limiting factor for performance. Using threads is automatic if you start Julia with multiple threads.

An experimental thread-parallel backend for matrices and linear algebra can be enabled by setting `backend=:csr` in the call to [`setup_reservoir_model`](@ref). This backend provides additional features such as a parallel zero-overlap ILU(0) implementation and parallel apply for AMG, but these features are still work in progress.

Starting Julia with multiple threads (for example `julia --project. --threads=4`) will allow `JutulDarcy` to make use of threads to speed up calculations

- The default behavior is to only speed up assembly of equations
- The linear solver is often the most expensive part -- as mentioned above, parts can be parallelized by choosing `csr` backend when setting up the model
- Running with a parallel preconditioner can lead to higher iteration counts since the ILU(0) preconditioner changes in parallel
- Heavy compositional models benefit a lot from using threads

Threads are easy to use and can give a bit of benefit for large models.

### Mixed-mode parallelism

You can mix the two approaches: Adding multiple threads to each MPI process can use threads to speed up assembly and property evaluations.

### Tips for parallel runs

A few hints when you are looking at performance:

- Reservoir simulations are memory bound, cannot expect that 10 threads = 10x performance
- CPUs can often boost single-core performance when resources are available
- MPI in JutulDarcy is less tested than single-process simulations, but is natural for larger models
- There is always some cost to parallelism: If running a large ensemble with limited compute, many serial runs handled by Julia's task system is usually a better option
- Adding the maximum number of processes does not always give the best performance. Typically you want at least 10 000 cells per process. Can be case dependent.

Example: 200k cell model on laptop: 1 process 235 s -> 4 processes 145s

## Solving with MPI in practice

There are a few adjustments needed before a script can be run in MPI.

### Setting up the environment

You will have to set up an environment with the following packages under Julia 1.9+:
`PartitionedArrays`, `MPI`, `JutulDarcy` and `HYPRE`. This is generally the best performing solver setup available, even if you are working in a shared memory environment.

### Writing the script

Write your script as usual. The following options must then be set:

- [`setup_reservoir_model`](@ref) should have the extra keyword argument `split_wells=true`. We also recommend `backend=:csr` for the best performance.
- [`simulate_reservoir`](@ref) or [`setup_reservoir_simulator`](@ref) should get the optional argument `mode = :mpi`

You must then run the file using the appropriate `mpiexec` as described in the MPI.jl documentation. Specialized functions will be called by `simulate_reservoir` when this is the case. We document them here, even if we recommend using the high level version of this interface:

```@docs
JutulDarcy.simulate_reservoir_parray
```

### Checklist for running in MPI

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

### Limitations of running in MPI

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
