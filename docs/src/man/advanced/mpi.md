# Multi-threading and MPI support

JutulDarcy can use threads by default, but advanced options can improve performance significantly for larger models.

## Overview of parallel support

There are two main ways of exploiting multiple cores in Jutul/JutulDarcy: Threads are automatically used for assembly and can be used for parts of the linear solve. If you require the best performance, you have to go to MPI where the linear solvers can use a parallel [BoomerAMG preconditioner](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) via [HYPRE.jl](https://github.com/fredrikekre/HYPRE.jl).

### MPI parallelization

MPI parallelizes all aspects of the solver using domain decomposition and allows a simulation to be divided between multiple nodes in e.g. a supercomputer. It is significantly more cumbersome to use than standard simulations as the program must be launched in MPI mode. This is typically a non-interactive process where you launch your MPI processes and once they complete the simulation the result is available on disk. The MPI parallel option uses a combination of MPI.jl, PartitionedArrays.jl and HYPRE.jl.

### Thread parallelization

JutulDarcy also supports threads. By defualt, this only parallelizes property evaluations and assembly of the linear system. For many problems, the linear solve is the limiting factor for performance.Using threads is automatic if you start Julia with multiple threads.

An experimental thread-parallel backend for matrices and linear algebra can be enabled by setting `backend=:csr` in the call to [`setup_reservoir_model`](@ref). This backend provides additional features such as a parallel zero-overlap ILU(0) implementation and parallel apply for AMG, but these features are still work in progress.

### Mixed-mode parallelism

You can mix the two approaches: Adding multiple threads to each MPI process can use threads to speed up assembly and property evaluations.

## Solving with MPI in practice

### Setting up the environment

Tou will have to set up an environment with the following packages under Julia 1.9+:
`PartitionedArrays`, `MPI`, `JutulDarcy` and `HYPRE`. This is generally the best performing solver setup available, even if you are working in a shared memory environment.

### Writing the script

Write your script as usual. The following options must then be set:

- [`setup_reservoir_model`](@ref) should have the extra keyword argument `split_wells=true`. We also recommend `backend=:csr` for the best performance.
- [`simulate_reservoir`](@ref) or [`setup_reservoir_simulator`](@ref) should get the optional argument `mode = :mpi`

You must then run the file using the approprioate `mpiexec` as described in the MPI.jl documentation. Specialized functions will be called by `simulate_reservoir` when this is the case. We document them here, even if we recommend using the high level version of this interface:

```@docs
setup_reservoir_simulator_parray
simulate_reservoir_parray
```

### Limitations for running in MPI

!!! note
    You should be familiar with the MPI programming model to use this feature. See [MPI.jl](https://juliaparallel.org/MPI.jl/stable/) for more details, and how MPI is handled in Julia specifically.

!!! note
    MPI consolidates results by writing files to disk. Unless you have a plan to work with the distributed states in-memory returned by the `simulate!` call, it is best to specify a `output_path` optional argument to [`setup_reservoir_simulator`](@ref). After the simulation, that folder will contain output just as if you had run the case in serial.
