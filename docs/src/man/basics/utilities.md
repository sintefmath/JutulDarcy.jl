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

## Adjoints and gradients

```@docs
JutulDarcy.reservoir_sensitivities
JutulDarcy.well_mismatch
Jutul.LBFGS.unit_box_bfgs
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
