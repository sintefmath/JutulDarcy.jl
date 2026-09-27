# GPU, multi-threading and MPI support

JutulDarcy uses a conservative configuration by default that is efficient for serial simulations, but you can employ different types of parallelism for speeding up larger models.

## Overview of parallel support

There are four main ways of running Jutul/JutulDarcy in parallel:

- Kernel-based threading. This parallelizes all aspects of the simulator and uses an alternative linear solver suitable for massive parallelization. **This is the preferred way of speeding up simulations on both CPU*. This can be enabled by passing `mode=:ka` to [`simulate_reservoir`](@ref) or [`setup_reservoir_simulator`](@ref).
- GPU execution. Massively parallelizes all parts of the code by running it on a GPU. See []
- MPI support. Run the same case across multiple compute nodes that do not share memory.
- Standard threading. This is automatically enabled by launching julia with threads (`julia --threads=N` where `N` is the number of threads). This speeds up many operations, but the default linear solver is only partially parallel. This is a conservative parallelization that works without changes to your simulation setup. You should likely try the kernel-based threading if you have large cases that can leverage many threads.

### MPI parallelization

MPI parallelizes all aspects of the solver using domain decomposition and allows a simulation to be divided between multiple nodes in e.g. a supercomputer. MPI is required to get the best parallel performance from the [BoomerAMG preconditioner](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) via [HYPRE.jl](https://github.com/fredrikekre/HYPRE.jl). It is significantly more cumbersome to use than standard simulations as the program must be launched in MPI mode. This is typically a non-interactive process where you launch your MPI processes and once they complete the simulation the result is available on disk. The MPI parallel option uses a combination of MPI.jl, PartitionedArrays.jl and HYPRE.jl.

For more details, see the [MPI support](@ref) page.

### Kernel parallelization (GPU/CPU)

This parallelizes all aspects of the simulator and is easy to use. We will cover the CPU case here.You can launch a kernel simulation by setting the mode:

```julia
simulate_reservoir(case, mode = :ka) # or mode = :ka_cpu
```

The GPU case is similar and covered in [GPU support](@ref) for CUDA/AMDGPU as the backends.

### Standard thread parallelization

JutulDarcy also supports threads. By default, this only parallelizes property evaluations and assembly of the linear system. For many problems, the linear solve is the limiting factor for performance. Using threads is automatic if you start Julia with multiple threads.

Starting Julia with multiple threads (for example `julia --project. --threads=4`) will allow `JutulDarcy` to make use of threads to speed up calculations

- The default behavior is to only speed up assembly of equations
- The linear solver is often the most expensive part -- as mentioned above, parts can be parallelized by choosing `csr` backend when setting up the model
- Running with a parallel preconditioner can lead to higher iteration counts since the ILU(0) preconditioner changes in parallel
- Heavy compositional models benefit a lot from using threads

Threads are easy to use and can give a bit of benefit for large models. You can also call `Jutul.set_hypre_threads(N)` where `N` is the desired number of threads to allow the hypre AMG solver to use threads as well. If you call `Jutul.set_hypre_threads()` it will default to the number of threads in the Julia session (i.e. whatever was passed to the `--threads` argument).

### Mixed-mode parallelism

You can mix MPI and threaded approaches approaches: Adding multiple threads to each MPI process can use threads to speed up assembly and property evaluations.

### Tips for parallel runs

A few hints when you are looking at performance:

- Reservoir simulations are memory bound, cannot expect that 10 threads = 10x performance
- CPUs can often boost single-core performance when resources are available
- MPI in JutulDarcy is less tested than single-process simulations, but is natural for larger models
- There is always some cost to parallelism: If running a large ensemble with limited compute, many serial runs handled by Julia's task system is usually a better option
- Adding the maximum number of processes does not always give the best performance. Typically you want at least 10 000 cells per process. Can be case dependent.

Example: 200k cell model on laptop: 1 process 235 s -> 4 processes 145s
