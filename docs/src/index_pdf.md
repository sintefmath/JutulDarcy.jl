# JutulDarcy.jl

Re-thinking reservoir simulation in Julia

High-performance porous media and reservoir simulator based on automatic differentiation

## What is this?

JutulDarcy.jl is a general high-performance purpose porous media simulator toolbox (CO2 sequestration, gas/H2 storage, oil/gas fields) written in [Julia](https://julialang.org/) based on [Jutul.jl](https://github.com/sintefmath/Jutul.jl), developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

A few highlights:

- Immiscible, black-oil, compositional, CO2-brine and geothermal systems
- Fully differentiable through adjoint method (history matching of parameters, optimization of well controls)
- High performance, with optional support for compiling MPI parallel binaries
- Consistent discretizations
- Industry standard input formats - or make your own model as a script
- 3D visualization and tools for post-processing of simulation results

## Key Features

### Physical systems
Immiscible, compositional, geothermal and black-oil flow

### Differentiability
Compute sensitivities of parameters with high-performance adjoint method

### High performance on CPU & GPU
Fast execution with support for MPI, CUDA and thread parallelism

## Quick start guide

[Getting started](@ref) is the main setup guide that includes the basics of installing Julia and creating a Julia environment for `JutulDarcy.jl`, written for users who may not already be familiar with Julia package management.

If you want to get started quickly: Install [Julia](https://julialang.org/) and add the following packages together
with a Makie backend library to your environment of choice using Julia's package manager `Pkg`:

```julia
using Pkg
Pkg.add("GLMakie")    # Plotting
Pkg.add("Jutul")      # Base package
Pkg.add("JutulDarcy") # Reservoir simulator
```

To verify that everything is working, we have a minimal example that runs an industry standard input file and produces interactive plots. Note that interactive plotting requires `GLMakie`, which may not work if you are running Julia over SSH.

### Python bindings

Alternatively, the code has a [Python package](https://github.com/sintefmath/PyJutulDarcy) that can be installed using `pip`:
```sh
pip install jutuldarcy
```

### Examples

The examples are published in the documentation. For a list of examples categorized by tags, see the Example overview page.

To get access to all the examples as code, you can generate a folder that contains the examples locally, you can run the following code to create a folder `jutuldarcy_examples` in your current working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples()
```

These examples can then be run using `include("jutuldarcy_examples/example_name.jl")` or opened in an editor to be run line by line.

## Citing JutulDarcy

If you use JutulDarcy for a scientific publication, please cite [the main paper](https://doi.org/10.3997/2214-4609.202437111) the following way:

> O. Møyner, (2024). JutulDarcy.jl - a Fully Differentiable High-Performance Reservoir Simulator Based on Automatic Differentiation. Computational Geosciences (2025), [Open Access: https://doi.org/10.1007/s10596-025-10366-6](https://doi.org/10.1007/s10596-025-10366-6)

BibTeX:

```bibtex
@article{jutuldarcy,
  title={JutulDarcy.jl - a fully differentiable high-performance reservoir simulator based on automatic differentiation},
  author={M{\o}yner, Olav},
  journal={Computational Geosciences},
  volume={29},
  number={30},
  year={2025},
  publisher={Springer}
}
```

## A few of the packages used by JutulDarcy

JutulDarcy.jl builds upon many of the excellent packages in the Julia ecosystem. Here are a few of them, and what they are used for:

- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) implements the Dual number class used throughout the code
- [SparsityTracing.jl](https://github.com/PALEOtoolkit/SparsityTracing.jl/) provides sparsity detection inside Jutul
- [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) provides the iterative linear solvers
- [ILUZero.jl](https://github.com/mcovalt/ILUZero.jl/blob/master/src/ILUZero.jl) for ILU(0) preconditioners
- [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) for AMG preconditioners
- [HYPRE.jl](https://github.com/fredrikekre/HYPRE.jl) for robust AMG preconditioners with MPI support
- [PartitionedArrays.jl](https://github.com/fverdugo/PartitionedArrays.jl) for MPI assembly and linear solve
- [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) for CUDA-GPU support
- [AMGX.jl](https://github.com/JuliaGPU/AMGX.jl) for AMG on CUDA GPUs
- [Tullio.jl](https://github.com/mcabbott/Tullio.jl) for automatically optimized loops and [Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl) for lightweight threads
- [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl) and [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) gives nice output to terminal.
- [Makie.jl](https://makie.juliaplots.org/) is used for the visualization features
- [MultiComponentFlash.jl](https://github.com/moyner/MultiComponentFlash.jl) provides many of the compositional features

...and many more directly, and indirectly - see the Project.toml and Manifest files for a full list!
