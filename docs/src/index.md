````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: JutulDarcy
  text:
  tagline: Reservoir simulation and porous media flow in Julia
  image:
    src: assets/logo.png
    alt: JutulDarcy
  actions:
    - theme: brand
      text: Getting started
      link: /man/intro
    - theme: alt
      text: View on Github
      link: https://github.com/sintefmath/JutulDarcy.jl
    - theme: alt
      text: API Reference
      link: /api
---
````

# Welcome to the JutulDarcy.jl documentation

JutulDarcy.jl: Darcy-scale and subsurface flow (CO2 sequestration, gas/H2 storage, oil/gas fields) using [Jutul.jl](https://github.com/sintefmath/Jutul.jl) developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Installation

Install [Julia](https://julialang.org/) and add the following packages together
with a Makie backend library to your environment of choice using Julia's package manager `Pkg`:

```julia
using Pkg
Pkg.add("GLMakie")    # Plotting
Pkg.add("Jutul")      # Base package
Pkg.add("JutulDarcy") # Reservoir simulator
```

You can then run any of the examples in the [`examples`](https://github.com/sintefmath/JutulDarcy.jl/tree/main/examples) directory by including them. The examples are sorted by complexity. We suggest you start with [Intro: Gravity segregation example](Intro: Gravity segregation example).

To generate a folder that contains the examples locally, you can run the following code to create a folder `jutuldarcy_examples` in your current working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples()
```

If you want a quick start, we also have a minimal example that runs a industry standard input file and produces interactive plots:

::: details Show me the code

```@example
using JutulDarcy, GLMakie
spe9_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE9")
file_path = joinpath(spe9_dir, "SPE9.DATA")
case = setup_case_from_data_file(file_path)
result = simulate_reservoir(case)
plot_reservoir_simulation_result(case.model, result)
```
:::

Note that interactive plotting requires `GLMakie`, which may not work if you are running Julia over SSH.

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
