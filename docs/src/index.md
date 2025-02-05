````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: JutulDarcy
  text: Re-thinking reservoir simulation in Julia
  image:
    src: logo.png
    alt: JutulDarcy
  tagline: High-performance porous media and reservoir simulator based on automatic differentiation
  actions:
    - theme: brand
      text: Getting started
      link: /man/intro
    - theme: alt
      text: View on Github
      link: https://github.com/sintefmath/JutulDarcy.jl
    - theme: alt
      text: Run .DATA file
      link: /examples/data_input_file

features:
  - icon: 💥
    title: Physical systems
    details: Immiscible, compositional, geothermal and black-oil flow
    link: /man/basics/systems

  - icon: ⚙️
    title: Differentiability
    details: Compute sensitivities of parameters with high-performance adjoint method
    link: /examples/intro_sensitivities

  - icon: 🏃
    title: High performance on CPU & GPU
    details: Fast execution with support for MPI, CUDA and thread parallelism
    link: /man/advanced/mpi
---
````

## What is this?

JutulDarcy.jl is a general high-performance purpose porous media simulator toolbox (CO2 sequestration, gas/H2 storage, oil/gas fields) written in [Julia](https://julialang.org/) based on [Jutul.jl](https://github.com/sintefmath/Jutul.jl), developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

A few highlights:

- Immiscible, black-oil, compositional, CO2-brine and geothermal systems
- Fully differentiable through adjoint method
- High performance, with optional support for compiling MPI parallel binaries
- Consistent discretizations
- Industry standard input formats - or make your own model as a script
- 3D visualization and tools for post-processing of simulation results

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

To verify that everything is working, we have a minimal example that runs a industry standard input file and produces interactive plots. Note that interactive plotting requires `GLMakie`, which may not work if you are running Julia over SSH.

```@example
using JutulDarcy, GLMakie
spe9_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE9")
file_path = joinpath(spe9_dir, "SPE9.DATA")
case = setup_case_from_data_file(file_path)
result = simulate_reservoir(case)
plot_reservoir_simulation_result(case.model, result)
```

To get access to all the examples, you can generate a folder that contains the examples locally, you can run the following code to create a folder `jutuldarcy_examples` in your current working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples()
```

These examples can then be run using `include("jutuldarcy_examples/example_name.jl")` or opened in an editor to be run line by line. Alternatively, you can download all examples as [Jupyter Notebooks](https://github.com/sintefmath/JutulDarcy.jl/tree/gh-pages/dev/final_site/notebooks).

## Citing JutulDarcy

If you use JutulDarcy for a scientific publication, please cite [the main paper](https://doi.org/10.3997/2214-4609.202437111) the following way:

> O. Møyner, (2024). JutulDarcy.jl - a Fully Differentiable High-Performance Reservoir Simulator Based on Automatic Differentiation. ECMOR 2024, [https://doi.org/10.3997/2214-4609.202437111](https://doi.org/10.3997/2214-4609.202437111)

::: details Show BibTeX

```bibtex
@article{jutuldarcy_ecmor_2024,
   author = "M{\o}yner, O.",
   title = "JutulDarcy.jl - a Fully Differentiable High-Performance Reservoir Simulator Based on Automatic Differentiation", 
   year = "2024",
   volume = "2024",
   number = "1",
   pages = "1-9",
   doi = "https://doi.org/10.3997/2214-4609.202437111",
   publisher = "European Association of Geoscientists \& Engineers",
   issn = "2214-4609",
}
```

:::
