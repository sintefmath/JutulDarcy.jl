# Parallel solution using MPI.jl, PartitionedArrays.jl and HYPRE.jl

JutulDarcy can use threads by default, but advanced options can improve performance significantly for larger models.

## Overview of parallel support

There are two main ways of exploiting multiple cores in Jutul/JutulDarcy: Threads are automatically used for assembly and can be used for parts of the linear solve. If you require the best performance, you have to go to MPI where the linear solvers can use a parallel [BoomerAMG preconditioner](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) via [HYPRE.jl](https://github.com/fredrikekre/HYPRE.jl).

### Shared memory

An experimental thread-parallel backend for matrices and linear algebra can be enabled by setting `backend=:csr` in the call to [`setup_reservoir_model`](@ref). This backend provides additional features such as a parallel zero-overlap ILU(0) implementation and parallel apply for AMG, but these features are still work in progress.

#### MPI support for distributed memory

It is possible to run cases using MPI. You will have to set up an environment with the following packages under Julia 1.9+:
`PartitionedArrays`, `MPI`, `JutulDarcy` and `HYPRE`. This is generally the best performing solver setup available, even if you are working in a shared memory environment.

### Why use MPI over threads?

### When to use threads instead of MPI?

### Mixed-mode parallelism

## Setting up environment for parallel solve

## Running the simulation

Write your script as usual, but in your call to `setup_reservoir_simulator`, pass the optional argument `mode = :mpi`. You must then run the file using `mpiexec`.

```@docs
setup_reservoir_simulator_parray
simulate_reservoir_parray
```

### MPI executable notes

### Limitations for running in MPI

!!! note
    You should be familiar with the MPI programming model to use this feature. See [MPI.jl](https://juliaparallel.org/MPI.jl/stable/) and [MPIClusterManagers.jl](https://github.com/JuliaParallel/MPIClusterManagers.jl) for more details, and how MPI is handled in Julia specifically.

!!! note
    MPI consolidates results by writing files to disk. Unless you have a plan to work with the distributed states in-memory returned by the `simulate!` call, it is best to specify a `output_path` optional argument to `setup_reservoir_simulator`. After the simulation, that folder will contain output just as if you had run the case in serial.
