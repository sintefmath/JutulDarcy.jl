[![DOI](https://zenodo.org/badge/477727603.svg)](https://zenodo.org/badge/latestdoi/477727603)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sintefmath.github.io/JutulDarcy.jl/dev/)
[![Build Status](https://github.com/sintefmath/JutulDarcy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sintefmath/JutulDarcy.jl/actions/workflows/CI.yml?query=branch%3Amain)


[![Jutul Darcy logo](https://github.com/sintefmath/JutulDarcy.jl/raw/main/docs/src/assets/logo_wide.png)](https://sintefmath.github.io/JutulDarcy.jl/dev/)

# Reservoir simulation in Julia
JutulDarcy.jl: Darcy-scale and subsurface flow (CO2 sequestration, gas/H2 storage, oil/gas fields) using [Jutul.jl](https://github.com/sintefmath/Jutul.jl) developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Key features
JutulDarcy is a general purpose porous media simulator with high performance written in Julia. It is fully differentiable with respect to forces and discretization parameters.

### Physical systems
- Immiscible multi-phase flow
- Black-oil type models with support for both dissolved vapor (Rs) and vaporized liquid (Rv)
- Equation-of-state compositional flow with up to three phases and any number of components

All solvers can incorporate general multisegment wells with rigorous mass balance, friction pressure loss, complex well limits and time-dependent controls.

### Architecture

- Written in pure Julia, with automatic differentiation and dynamic sparsity detection
- Support for sensitivities with respect to any model parameters using the adjoint method
- High performance assembly and linear solvers, with support for two-stage CPR BILU(0)-CPR Krylov solvers
- MPI support with domain decomposition and BoomerAMG-CPR solver with automatic METIS partitioning

### Input/output and workflow tools
- Unstructured grids and complex cases input from [the Matlab Reservoir Simulation Toolbox (MRST)](https://www.mrst.no) using the `jutul` module
- Experimental support for reading in and running .DATA files with corner point grids
- 3D visualization of grids and wells by loading a [Makie.jl](https://docs.makie.org/stable/) backend (requires Julia 1.9, `GLMakie` for interactivity)
- Interactive plotting of well curves

The compositional simulator has been matched against commercial offerings, AD-GPRS and MRST. The blackoil simulator has been validated on the standard SPE benchmarks (SPE1, SPE9, ...).

## Example run times on benchmarks
| Name      | Cells | Report steps | Preconditioner   | Time [s] |
|-----------|-------|--------------|------------------|----------|
| SPE1CASE2 | 300   | 120          | block-ILU(0)     | 0.30     |
| SPE9      | 9000  | 35           | block-ILU(0)     | 3.41     |
| Egg       | 18553 | 123          | CPR-block-ILU(0) | 8.60     |

Simulated with default `julia`, no threads. Runtime with single-threaded HYPRE on the full [SPE10 model 2 model](https://www.spe.org/web/csp/datasets/set02.htm) is about 20 minutes with default settings and strict mass-balance tolerances and less than 4 minutes for realization 1 of the [OLYMPUS optimization benchmark challenge](https://link.springer.com/article/10.1007/s10596-020-10003-4).

## A few of the packages used by Jutul and JutulDarcy
Jutul builds upon many of the excellent packages in the Julia ecosystem. Here are a few of them, and what they are used for:
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) implements the Dual number class used throughout the code
- [SparsityTracing.jl](https://github.com/PALEOtoolkit/SparsityTracing.jl/) provides sparsity detection inside Jutul
- [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) provides the iterative linear solvers
- [ILUZero.jl](https://github.com/mcovalt/ILUZero.jl/blob/master/src/ILUZero.jl) for ILU(0) preconditioners
- [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) for AMG preconditioners
- [HYPRE.jl](https://github.com/fredrikekre/HYPRE.jl) for robust AMG preconditioners with MPI support
- [PartitionedArrays.jl](https://github.com/fverdugo/PartitionedArrays.jl) for MPI assembly and linear solve
- [Tullio.jl](https://github.com/mcabbott/Tullio.jl) for automatically optimized loops and [Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl) for lightweight threads
- [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl) and [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) gives nice output to terminal 
- [Makie.jl](https://makie.juliaplots.org/) is used for the visualization features
- [MultiComponentFlash.jl](https://github.com/moyner/MultiComponentFlash.jl) provides many of the compositional features

...and many more, both directly in the Project.toml file and indirectly!

## Getting started
Install [Julia](https://julialang.org/) and add the package to your environment of choice:
```julia
using Pkg
Pkg.add("CairoMakie")
Pkg.add("Jutul")
Pkg.add("JutulDarcy")
```
You can then run any of the examples in the `examples` directory by including them.

# Additional examples and further reading
The [documentation](https://sintefmath.github.io/JutulDarcy.jl/dev/) is still work in progress, but contains a fair bit of useful information. In addition, see the [examples](examples/) folder for more information. Some functionality is also demonstrated in the [tests](https://github.com/sintefmath/JutulDarcy.jl/tree/main/test).

*Internals and undocumented functions are subject to change at this time. However, the main interface for the reservoir simulator itself seen in the examples should be fairly stable.*
