# Standalone reservoir simulator

Scripts are interactive and useful for doing setup, simulation and post-processing in one file, but sometimes you want to run a big model unmodified from an input file:

- As an alternative to a pure Julia workflow, `JutulDarcy.jl` can be compiled into a standalone reservoir simulator
- This makes MPI simulations more ergonomic
- Compiling the code saves time when running multiple simulations
- The resulting executable is a standard command-line program - no Julia experience needed
- Output is given in the same format as regular simulations, can load data by restarting a simulation from the same `output_path`

This workflow uses [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl). For more details and an example build file with keyword arguments, see [the JutulDarcyApps.jl repository](https://github.com/sintefmath/JutulDarcyApps.jl/tree/master/mpi_simulator).
