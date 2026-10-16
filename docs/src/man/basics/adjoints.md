# Adjoints and gradients

JutulDarcy models are fully differentiable, and can produce gradients with
respect to any input parameter through the adjoint method. For more details, we
refer you to the examples. It is recommended to go through the [Introduction to
history matching with gradients](@ref) example to get an overview of the
functionality.

## Calculation of gradients

```@docs
JutulDarcy.reservoir_sensitivities
JutulDarcy.parameters_gradient_reservoir
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
series). This amounts to first setting up an observation that computes (typically) sum of squares mismatch with weighting for the wells you want to match, and then parametrizing the model using either a custom setup function or `setup_reservoir_dict_optimization`.

```@docs
JutulDarcy.HistoryMatching.history_match_objective
JutulDarcy.HistoryMatching.match_well!
JutulDarcy.HistoryMatching.match_injectors!
JutulDarcy.HistoryMatching.match_producers!
```

### Advise for history matching

If you are setting up a model for matching, we make the following recommendations for best practice. The tips here are general and can be adapted if you have a geological concept in mind.

- Group layers in the model for permeability and porosity by a coarse lateral position and layers in the model. This is either done by writing a custom setup function, or by using the `lumping` keyword argument to `free_optimization_parameter`. Scale permeability by log (`scaler = :log` in the same call).
- Set reasonable bounds for the permeability and porosity.
- Set vertical/horizontal permeability ratio to vary.
- Set `optimizer = :lbfgs_qp` when running `optimize_reservoir`.
- Let the producers operate under a total liquid rate, and then match the oil and gas rates.
- Balance the weights so that the contributions between the rates and the bhp are approximately proportional to the initial error (can be plotted using `JutulDarcy.plot_mismatch`).
- Let the optimizer run for many iterations if you can afford the runtime (e.g. 50 or so iterations).
- Use accelerated modes (pass e.g. `mode = :ka` to the `optimize_reservoir` call) if you have many cells and threads enabled. See also [GPU support](@ref).

## Utilities

```@docs
JutulDarcy.setup_reservoir_dict_optimization
JutulDarcy.well_mismatch
JutulDarcy.compute_well_qoi
JutulDarcy.npv_objective
```
