# Adjoints and gradients

JutulDarcy models are fully differentiable, and can produce gradients with
respect to any input parameter through the adjoint method. For more details, we
refer you to the examples. It is recommended to go through the [Introduction to
history matching with gradients](@ref) example to get an overview of the
functionality.

## Calculation of gradients

```@docs
JutulDarcy.reservoir_sensitivities
```

## Optimization interface

We use the functions from `Jutul` to free parameters before optimization.

```@docs
Jutul.free_optimization_parameter!
Jutul.add_optimization_multiplier!
optimize_reservoir
Jutul.DictOptimization.DictParameters
```

## Functionality from Jutul.jl

This optimization interface is a wrapper around the Jutul function `optimize` with preselected options that are sensible for reservoir models. It can still be useful to look at the inner function to see additional supported arguments.

```@docs
Jutul.optimize
Jutul.LBFGS.unit_box_bfgs
```

## History matching module

The history matching module makes it easy to define objective functions from
typical reservoir observations (rates, pressures and fractions per well as time
series).

```@docs
JutulDarcy.HistoryMatching.history_match_objective
JutulDarcy.HistoryMatching.match_well!
JutulDarcy.HistoryMatching.match_injectors!
JutulDarcy.HistoryMatching.match_producers!
```

## Utilities

```@docs
JutulDarcy.setup_reservoir_dict_optimization
JutulDarcy.well_mismatch
JutulDarcy.compute_well_qoi
JutulDarcy.npv_objective
```
