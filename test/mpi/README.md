Run the MPI regression test from the repository root:

```sh
julia --project=test/mpi -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=test/mpi test/mpi/run.jl
```

The launcher starts two Julia processes with MPI.jl's matching `mpiexec` and a
shared temporary output directory. The test requires ghost rows on both ranks,
uses CSR matrices and CPR/HYPRE, and compares the consolidated reservoir and well
results with a serial simulation. CI runs it separately from the serial suite.
