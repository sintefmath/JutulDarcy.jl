# Frequently asked questions

Here are a few common questions and possible answers. You may also want to have a look at the [GitHub issues](https://github.com/sintefmath/JutulDarcy.jl/issues) and [the GitHub discussions page](https://github.com/sintefmath/JutulDarcy.jl/discussions).

## Input and output

### What input formats can JutulDarcy.jl use?

1. DATA files (used by Eclipse, OPM Flow, tNavigator, Echelon and others) provided that the grid is given either as a corner-point GRDECL file or in TOPS format. As with most reservoir simulators, not all features of the original format are supported, but the code will let you know when unsupported features are encountered.
1. AFI-files (used by Intersect), provided that the definition of the grid and grid quantities is either given inline or embedded as RESQML. This format is more experimental than DATA files as less samples are publicly available.
1. Cases written out from [MRST](https://www.mrst.no/) through the `jutul` module.
1. Cases written entirely in Julia using the basic `Jutul` and `JutulDarcy` data structures, as seen in the examples of the module, with options to use other packages for e.g. mesh generation or data import.

### What output formats does JutulDarcy.jl have?

The simulator outputs results into standard Julia data structures (e.g. Vectors and Dicts) that can easily be written out using other Julia packages, for example in CSV format. We do not currently support binary formats output by commercial simulators.

Simulation results are written to disk using [JLD2](https://github.com/JuliaIO/JLD2.jl), a subset of [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) commonly used in Julia for storing objects to disk.

### What units does JutulDarcy.jl use?

JutulDarcy nominally uses SI units, but the code contains utilities to convert from typical unit sets to the internal SI format. See [Unit support](@ref) for more details and a tutorial. Wrong units (caused by e.g. directly inserting field units like millidarcy or standard barrels without conversion) can lead to large differences between simulators, or make problems very difficult to simulate.

### How do I restart an interrupted simulation?

JutulDarcy keeps everything in memory by default. This is not practical for larger models. If the argument `output_path` is set to a directory, JutulDarcy writes to the `JLD2` format (variant of HDF5).

```julia
# Note: set ENV["JUTUL_OUTPUT_PATH"] in your startup.jl first!
pth = jutul_output_path("My_test_case")
simulate_reservoir(case, output_path = pth)
```

If an output path is set, you can restart simulations:

```julia
# Restart from the last successfully solved step, or return output if everything is simulated
ws, states = simulate_reservoir(case, output_path = pth, restart = true)
# Start from the beginning, overwriting files if already present
ws, states = simulate_reservoir(case, output_path = pth, restart = false)
# Restart from step 10 and throw error if step 9 is not already stored on disk.
ws, states = simulate_reservoir(case, output_path = pth, restart = 10)
```

You can restart the simulation with different options for timestepping or tolerances.

### How do I decide where output is stored?

`Jutul.jl` contains a system for managing output folders. It is highly recommended that you amend your `startup.jl` file to include `ENV["JUTUL_OUTPUT_PATH"]` that points to where you want output to be stored. For example, on Windows usage of the output path mechanism may look something like this:

```julia
julia> ENV["JUTUL_OUTPUT_PATH"]
"D:/jutul_output/"

julia> jutul_output_path() # Randomly generated file name
"D:/jutul_output/jutul/jl_DwpAvQTiLo"

julia> jutul_output_path("mycase")
"D:/jutul_output/jutul/mycase"

julia> jutul_output_path("mycase", subfolder = "ensemble_name")
"D:/jutul_output/ensemble_name/mycase"

julia> jutul_output_path("mycase", subfolder = missing)
"D:/jutul_output/mycase"
```

Or equivalent on a Linux system:

```julia
julia> ENV["JUTUL_OUTPUT_PATH"]
"/home/username/jutul_output/"

julia> jutul_output_path() # Randomly generated file name
"/home/username/jutul_output/jutul/jl_DwpAvQTiLo"

julia> jutul_output_path("mycase")
"/home/username/jutul_output/jutul/mycase"

julia> jutul_output_path("mycase", subfolder = "ensemble_name")
"/home/username/jutul_output/ensemble_name/mycase"

julia> jutul_output_path("mycase", subfolder = missing)
"/home/username/jutul_output/mycase"
```

You can also just specify a full path and keep track of output folders yourself, but using the `jutul_output_path` mechanism will make it easier to write a script that can be run on another computer with different folder structure.

### How do you get out more output from a simulation?

The default outputs per cell are primary variables and total masses:

```julia
reservoir_model(model).output_variables
3-element Vector{Symbol}:
 :Pressure
 :Saturations
 :TotalMasses
```

You can push variables to this list, or ask the code to output all variables:

```julia
model2, = setup_reservoir_model(domain, sys, extra_outputs = true);
reservoir_model(model2).output_variables
7-element Vector{Symbol}:
 :Pressure
 :Saturations
 :TotalMasses
 :PhaseMassDensities
 :RelativePermeabilities
 :PhaseMobilities
 :PhaseMassMobilities
```

You can also pass `extra_outputs = [:PhaseMobilities]` as a keyword argument to `setup_reservoir_model` to make the resulting model output specific variables.

### What is the unit and sign convention for well rates?

Well results are given in strict SI, which means that rates are generally given in ``m^3/s``. Rates are positive for injection (mass entering the reservoir domain) and negative for production (leaving the reservoir domain).

## Plotting

### What is required for visualization?

We use the wonderful [`Makie.jl`](https://docs.makie.org/) for both 2D and 3D plots. Generally `CairoMakie` is supported for non-interactive plotting and `GLMakie` is used for interactive plotting (especially 3D). The latter requires a working graphics context, which is not directly available when the code is run over for example SSH or on a server.

For more details on the backends, see the[`Makie.jl` docs](https://docs.makie.org/stable/explanations/backends/backends)

## Miscellaneous

### Can you add feature X or format Y?

If you have a feature you'd like to have supported, please file an
[issue](https://github.com/sintefmath/JutulDarcy.jl/issues) with details on the format.
`JutulDarcy` is developed primarily through contract research, so features are
added as needed for ongoing projects where the simulator is in use. Posting an
issue, especially if you have a clear reference to how something should be
implemented is still very useful. It is also possible to fund the development
for a specific feature, or to implement the feature yourself by asking for
pointers on how to get started.

### Who develops JutulDarcy.jl?

The module is developed and maintained by the [Applied Computational Sciences group](https://www.sintef.no/en/digital/departments-new/department-of-mathematics-and-cybernetics/research-group-applied-computational-science/) at [SINTEF Digital](https://www.sintef.no/sintef-digital/). SINTEF is one of Europe's largest independent research organizations and is organized as a not-for-profit institute. [Olav Møyner](https://www.sintef.no/en/all-employees/employee/olav.moyner/) is the primary maintainer.

### Why write a new reservoir simulation code in Julia?

We believe that reservoir simulation should be a *library* and not necessarily an application by itself. The future of porous media simulation is deeply integrated into other workflows, and not as an application that simply writes files to disk. As a part of experimentation in differentiable and flexible solvers using automatic differentiation that started with [MRST](https://www.mrst.no), Julia was the natural next step.

### What is the license of JutulDarcy.jl?

The code uses the [MIT license](https://en.wikipedia.org/wiki/MIT_License), a permissive license that requires attribution, but does not place limitations on commercial use or closed-source integration.

The code uses a number of dependencies that can have other licenses and we make no guarantees that the entirety of the code made available by adding `JutulDarcy.jl` to a given Julia environment is all MIT licensed.
