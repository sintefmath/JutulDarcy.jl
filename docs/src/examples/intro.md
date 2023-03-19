## Getting started
Install packages:
```julia
using Pkg
Pkg.add("Jutul")
Pkg.add("JutulDarcy")
```

If you want the plotting used in the examples, you need this:
```julia
Pkg.add("CairoMakie")
```

In addition, there is experimental 3D visualization and well plotting found in an unregistered package. At some point in time this will be registered or folded into the main package as a conditional dependency.

```julia
Pkg.add(url="https://github.com/sintefmath/JutulViz.jl")
```
