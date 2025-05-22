
# Standalone reservoir simulator {#Standalone-reservoir-simulator}

Scripts are interactive and useful for doing setup, simulation and post-processing in one file, but sometimes you want to run a big model unmodified from an input file:
- As an alternative to a pure Julia workflow, `JutulDarcy.jl` can be compiled into a standalone reservoir simulator
  
- This makes MPI simulations more ergonomic
  
- Compiling the code saves time when running multiple simulations
  
- The resulting executable is a standard command-line program - no Julia experience needed
  
- Output is given in the same format as regular simulations, can load data by restarting a simulation from the same `output_path`
  

This workflow uses [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl). For more details and an example build file with keyword arguments, see [the JutulDarcyApps.jl repository](https://github.com/sintefmath/JutulDarcyApps.jl/tree/master/mpi_simulator).

A few things to note:
- The simulator comes with a set of shared library files and will be ~500 mb
  
- Binaries will match platform (compiling under Linux gives you Linux binaries)
  
- The repository has a script that runs small &quot;representative&quot; models
  
- You can input small representative models in `precompile_jutul_darcy_mpi.jl` to make sure that compilation is avoided during simulation
  
- By default, the script uses the default Julia MPI binary. On a cluster, the build script may have to be modified to use the MPI type of the cluster using [MPITrampoline.jl](https://github.com/eschnett/MPItrampoline)
  

If you get it working on a complex MPI setup, feedback on your experience and PRs are very welcome.
