[![Paper](https://img.shields.io/badge/open_access_paper-computational_geosciences-blue.svg)](https://link.springer.com/article/10.1007/s10596-025-10366-6)
[![Zenodo](https://zenodo.org/badge/477727603.svg)](https://zenodo.org/badge/latestdoi/477727603)
[![Dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://sintefmath.github.io/JutulDarcy.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sintefmath.github.io/JutulDarcy.jl/dev/)
[![Build Status](https://github.com/sintefmath/JutulDarcy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sintefmath/JutulDarcy.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FJutulDarcy&query=total_requests&suffix=%2Fmonth&label=Downloads&color=green)](https://juliapkgstats.com/pkg/JutulDarcy)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FJutulDarcy&query=total_requests&&label=Total%20Downloads&color=green)](https://juliapkgstats.com/pkg/JutulDarcy)


[![Jutul Darcy logo](https://github.com/sintefmath/JutulDarcy.jl/raw/main/docs/src/assets/logo_wide.png)](https://sintefmath.github.io/JutulDarcy.jl/stable/)

> [!TIP]
> Visit the docs at https://sintefmath.github.io/JutulDarcy.jl/stable/

# Reservoir simulation in Julia

JutulDarcy.jl is a toolkit for building, running and interacting with reservoir simulation models written in Julia. The software covers Darcy-scale subsurface flow with wells (CO2 sequestration, geothermal reservoirs, gas/H2 storage, oil/gas fields) using the [Jutul.jl](https://github.com/sintefmath/Jutul.jl) framework for differentiable programming developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Key features

JutulDarcy is a general purpose porous media simulator toolkit with high performance written in Julia. It is fully differentiable with respect to forces and discretization parameters. The package contains utilities for building, plotting, running and analysing models. Extensive [validation examples](https://sintefmath.github.io/JutulDarcy.jl/dev/man/validation) provide comparison to other solvers.

### Physical systems

- Single-phase flow
- Immiscible multi-phase flow
- Black-oil type models with support for both dissolved vapor (Rs) and vaporized liquid (Rv)
- Equation-of-state compositional flow with up to three phases and any number of components (Peng-Robinson, Søreide-Whitson, Soave-Redlich Kwong, ...)
- Specialized CO2-brine model for CO2 storage
- Thermal effects
- Tracers and polymer EOR
- Geothermal systems (see also [Fimbul.jl](https://github.com/sintefmath/Fimbul.jl) for extended geothermal capabilities)

All solvers can incorporate general multisegment wells with rigorous mass balance, friction pressure loss, complex well limits and time-dependent controls.

### Architecture

- Written in pure Julia, with automatic differentiation and dynamic sparsity detection
- Support for sensitivities with respect to any model parameters using the adjoint method
- Support for GPU execution via KernelAbstractions - run in parallel on both NVIDIA and AMD cards.
- High performance assembly and linear solvers, with support for two-stage CPR BILU(0)-CPR Krylov solvers
- MPI support with domain decomposition and BoomerAMG-CPR solver with automatic METIS partitioning
- Support for consistent discretizations (AvgMPFA / NTPFA)

### Input/output and workflow tools

- Support for [reading in and running .DATA and .AFI files](https://sintefmath.github.io/JutulDarcy.jl/stable/man/basics/input_files#DATA-files-from-commercial-reservoir-modelling-software) with corner point grids, with support for non-conformal faults and inactive cells.
- Unstructured grids and complex cases input from [the Matlab Reservoir Simulation Toolbox (MRST)](https://www.mrst.no) using the `jutul` module
- 3D visualization of grids and wells by loading a [Makie.jl](https://docs.makie.org/stable/) backend (requires Julia 1.10, `GLMakie` for interactivity)
- Interactive plotting of well curves, field responses and other quantities of interest.
- Automatically coarsen models to reduce runtime.

## Example run times on benchmarks

### Small benchmarks

| Name      | Cells  | Report steps | Preconditioner   | Time [s] |
|-----------|--------|--------------|------------------|----------|
| SPE1CASE2 | 300    | 120          | block-ILU(0)     | 0.30     |
| SPE9      | 9000   | 35           | block-ILU(0)     | 3.41     |
| Egg       | 18553  | 123          | CPR-block-ILU(0) | 8.60     |
| Norne     | 44417  | 247          | CPR-block-ILU(0) | 228.8    |

 Numbers are given for serial CPU solves with `hypre` as the pressure solver in CPR.

### Parallel/GPU benchmarks

JutulDarcy supports parallelism by leveraging threads, MPI and GPU (CUDA/AMD).
The solvers are fully resident on GPU, allowing for much faster solves for
larger models. A few examples on consumer grade hardware (Ryzen 9 7950X with 15
threads, and a RTX4080 GPU):

| Name | Cells | Type | CPU Time [s] | GPU Time [s] |
|---|---:|---|---:|---:|
| [SPE10 model 2](https://www.spe.org/web/csp/datasets/set02.htm) | 1,082,465 | 2-phase immiscible | 179.0 | 11.8 |
| [OLYMPUS realization 1](https://link.springer.com/article/10.1007/s10596-020-10003-4) | 192,750 | 2-phase immiscible | 75.8 | 32.5 |
| [UHS Benchmark study](https://www.ite.tu-clausthal.de/en/research/subsurface-energy-and-gas-storage/uhs-benchmark-study) | 44,652 | 8-component EoS | 170.4 | 56.3 |

## Paper and citing

The main paper describing `JutulDarcy.jl` is [*JutulDarcy.jl - a Fully Differentiable High-Performance Reservoir Simulator Based on Automatic Differentiation*, available as open access from Springer](https://link.springer.com/article/10.1007/s10596-025-10366-6):

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

[DOI link to open access paper.](https://doi.org/10.1007/s10596-025-10366-6) If you use JutulDarcy in your work, please cite this paper.

[List of papers that make use of JutulDarcy.jl (add your own!)](https://sintefmath.github.io/JutulDarcy.jl/dev/extras/paper_list)

## Getting started

### Julia

Install [Julia](https://julialang.org/) and add the package to your environment of choice:

```julia
using Pkg
Pkg.add("GLMakie")
Pkg.add("Jutul")
Pkg.add("JutulDarcy")
```

You can then run any of the examples in the [`examples`](https://github.com/sintefmath/JutulDarcy.jl/tree/main/examples) directory by including them.

For a minimal example that runs an industry standard input file and produces interactive plots:

```julia
using JutulDarcy, GLMakie
spe9_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE9")
file_path = joinpath(spe9_dir, "SPE9.DATA")
case = setup_case_from_data_file(file_path)
result = simulate_reservoir(case)
plot_reservoir_simulation_result(case.model, result)
```

Note that interactive plotting requires `GLMakie`, which may not work if you are running Julia over SSH.

### Python

There is also a Python interface. You can install JutulDarcy using pip:

```bash
pip install jutuldarcy
```

See the [PyJutulDarcy documentation](https://sintefmath.github.io/PyJutulDarcy/) for setup, examples, and the Python API reference.

### Compiled executable

It is also possible to build an executable. See [the JutulDarcyApps.jl repository](https://github.com/sintefmath/JutulDarcyApps.jl/tree/master/mpi_simulator) for more information.

# Additional examples and further reading

The [documentation](https://sintefmath.github.io/JutulDarcy.jl/dev/) is the main source of information, [including published versions of the examples](https://sintefmath.github.io/JutulDarcy.jl/dev/examples/introduction/wells_intro). In addition, see the [examples](https://github.com/sintefmath/JutulDarcy.jl/tree/main/examples) folder for more information. Some functionality is also demonstrated in the [tests](https://github.com/sintefmath/JutulDarcy.jl/tree/main/test).

*Internals and undocumented functions are subject to change at this time. However, the main interface for the reservoir simulator itself is stable and follows semantic versioning for breakage.*
