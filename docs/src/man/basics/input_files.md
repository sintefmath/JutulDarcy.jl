
# Input formats

It is also possible to read cases that have been set up in MRST (see [`setup_case_from_mrst`](@ref) and [`simulate_mrst_case`](@ref)) or from .DATA files (see [`setup_case_from_data_file`](@ref) and [`simulate_data_file`](@ref))

## MAT-files from the Matlab Reservoir Simulation Toolbox (MRST)

### Simulation of .MAT files

```@docs
setup_case_from_mrst
simulate_mrst_case
```

### MRST-specific types

```@docs
Jutul.MRSTWrapMesh
```

## DATA-files from commercial reservoir modelling software

JutulDarcy can set up cases from Eclipse-type input files by making use of the [GeoEnergyIO.jl](https://github.com/sintefmath/GeoEnergyIO.jl) package for parsing. This package is a direct dependency of JutulDarcy and these cases can be simulated directly. If you want to parse the input files and possibly modify them in your Julia session before the case is simulated, we refer you to the [GeoEnergyIO Documentation](https://sintefmath.github.io/GeoEnergyIO.jl/dev/).

If you want to directly simulate a file from disk, you can sue the high level functions that automatically parse the files for you:

```@docs
simulate_data_file
setup_case_from_data_file
JutulDarcy.setup_case_from_parsed_data
JutulDarcy.convert_co2store_to_co2_brine
```

Reservoir simulator input files are highly complex and contain a great number of
different keywords. JutulDarcy and GeoEnergyIO has extensive support for this
format, but many keywords are missing or only partially supported. Sometimes
cases can be made to run by removing keywords that do not impact simulation
outcomes, or have very little impact. If you want a turnkey open source solution
for simulation reservoir models you should also have a look at [OPM
Flow](https://opm-project.org/).

If you can share your input file or the missing keywords [in the issues
section](https://github.com/sintefmath/JutulDarcy.jl/issues) it may be easier to
figure out if a case can be supported.
