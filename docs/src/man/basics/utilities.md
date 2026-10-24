# Utilities

This section describes various utilities that do not fit in other sections.

## CO2 and brine correlations

These functions are not exported, but can be found inside the `CO2Properties` submodule. The functions described here are a Julia port of the MRST module described in [salo_co2](@cite) They can be accessed by explicit import:

```julia
import JutulDarcy.CO2Properties: name_of_function
```

```@docs
JutulDarcy.CO2Properties.compute_co2_brine_props
JutulDarcy.CO2Properties.pvt_brine_RoweChou1970
JutulDarcy.CO2Properties.activity_co2_DS2003
JutulDarcy.CO2Properties.viscosity_brine_co2_mixture_IC2012
JutulDarcy.CO2Properties.pvt_brine_BatzleWang1992
JutulDarcy.CO2Properties.viscosity_co2_Fenghour1998
JutulDarcy.CO2Properties.pvt_co2_RedlichKwong1949
JutulDarcy.CO2Properties.viscosity_gas_mixture_Davidson1993
```

## Relative permeability functions

```@docs
JutulDarcy.brooks_corey_relperm
```

## CO2 inventory

```@docs
JutulDarcy.co2_inventory
JutulDarcy.plot_co2_inventory
```

## Geostatistical realizations

Spatially correlated porosity and permeability fields can be generated with
[GeoStats.jl](https://juliaearth.github.io/GeoStats.jl) once it is loaded
(`using GeoStats`). For meshes with a logical IJK structure (Cartesian and
corner-point meshes) the realization is by default sampled on the logical grid,
so that the correlation follows the layers of the model. Properties defined on
a separate Cartesian grid can be mapped to any mesh with volume-weighted
averaging, both when the grid is finer and when it is coarser than the mesh.

```julia
using JutulDarcy, GeoStats
domain = reservoir_domain(mesh)
realization = generate_perm_poro!(domain, seed = 1, horizontal_correlation_range = 200.0, vertical_correlation_range = 5.0)
# Map values from a fine Cartesian grid, e.g. a GeoStats realization
map_property_to_cells!(domain, grid, Dict(:porosity => poro, :permeability => perm))
```

```@docs
generate_perm_poro
generate_perm_poro!
map_property_to_cells
map_property_to_cells!
kozeny_carman_permeability
```

## API utilities

```@docs
reservoir_model
JutulDarcy.reservoir_storage
JutulDarcy.well_symbols
```

## Model reduction

```@docs
coarsen_reservoir_case
```

## Well outputs

```@docs
get_model_wells
well_output
full_well_outputs
```

## Non-neighboring connections

```@docs
setup_nnc_connections
```
