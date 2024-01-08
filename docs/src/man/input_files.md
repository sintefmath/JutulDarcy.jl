
# Reading input files

It is also possible to read cases that have been set up in MRST (see [`setup_case_from_mrst`](@ref) and [`simulate_mrst_case`](@ref)) or from .DATA files (see [`parse_data_file`](@ref) and [`simulate_data_file`](@ref))

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

### Parsers

```@docs
parse_data_file
parse_grdecl_file
```

### Simulation of .DATA files

```@docs
simulate_data_file
```
