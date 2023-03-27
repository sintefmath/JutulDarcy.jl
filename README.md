[![DOI](https://zenodo.org/badge/477727603.svg)](https://zenodo.org/badge/latestdoi/477727603)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sintefmath.github.io/JutulDarcy.jl/dev/)
[![Build Status](https://github.com/sintefmath/JutulDarcy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sintefmath/JutulDarcy.jl/actions/workflows/CI.yml?query=branch%3Amain)


![Jutul Darcy logo](docs/src/assets/logo_wide.png)

# Reservoir simulation in Julia
JutulDarcy.jl: Darcy-scale and subsurface flow (CO2 sequestration, gas/H2 storage, oil/gas fields) using [Jutul.jl](https://github.com/sintefmath/Jutul.jl) developed by the [Computational Geosciences group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/computational-geoscience/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Key features
- Written in pure Julia, with automatic differentiation and dynamic sparsity detection
- Support for sensitivities with respect to any model parameters using the adjoint method
- High performance assembly and linear solvers, with support for two-stage CPR BILU(0)-CPR Krylov solvers
- Equation-of-state compositional, immiscible and black oil flow is supported and validated against existing simulators
- Unstructured grids and complex cases input from [the Matlab Reservoir Simulation Toolbox (MRST)](https://www.mrst.no) using the `jutul` module.
- Support for general multisegment wells with rigorous mass balance, complex well limits and time-dependent controls
- 3D visualization of grids and wells in [JutulViz.jl](https://github.com/sintefmath/JutulViz.jl)
- Interactive plotting of well curves

The compositional simulator has been matched against commercial offerings, AD-GPRS and MRST. The blackoil simulator has been validated on the standard SPE benchmarks (SPE1, SPE9, ...).

## Example run times on benchmarks
| Name      | Cells | Report steps | Preconditioner   | Time [s] |
|-----------|-------|--------------|------------------|----------|
| SPE1CASE2 | 300   | 120          | block-ILU(0)     | 0.30     |
| SPE9      | 9000  | 35           | block-ILU(0)     | 3.41     |
| Egg       | 18553 | 123          | CPR-block-ILU(0) | 8.60     |

Simulated with `julia -O2`, no threads.

## A few of the packages used by Jutul and JutulDarcy
Jutul builds upon many of the excellent packages in the Julia ecosystem. Here are a few of them, and what they are used for:
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) implements the Dual number class used throughout the code
- [SparsityTracing.jl](https://github.com/PALEOtoolkit/SparsityTracing.jl/) provides sparsity detection inside Jutul
- [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) provides the iterative linear solvers
- [ILUZero.jl](https://github.com/mcovalt/ILUZero.jl/blob/master/src/ILUZero.jl) for ILU(0) preconditioners
- [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) for AMG preconditioners
- [Tullio.jl](https://github.com/mcabbott/Tullio.jl) for automatically optimized loops and [Polyester.jl]() for lightweight threads
- [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl) and [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) gives nice output to terminal 
- [Makie.jl](https://makie.juliaplots.org/) is used for the visualization features found in [JutulViz.jl](https://github.com/sintefmath/JutulViz.jl)
- [MultiComponentFlash.jl](https://github.com/moyner/MultiComponentFlash.jl) provides many of the compositional features

...and many more, both directly in the Project.toml file and indirectly!

## Getting started
Install [Julia](https://julialang.org/) and add the package to your environment of choice:
```julia
using Pkg
Pkg.add("Jutul")
Pkg.add("JutulDarcy")
Pkg.add("CairoMakie")
```
You can then run any of the examples in the `examples` directory by including them.

# Additional examples and further reading
As per Julia tradition, the documentation is fairly sparse at the moment. Please see the [examples](examples/) folder for more information. There are also examples in [JutulExamples](https://github.com/sintefmath/JutulExamples.jl) that may occasionally be a bit out of date.

*Internals and undocumented functions are subject to change at this time. However, the main interface for the reservoir simulator itself seen in the examples should be fairly stable.*
