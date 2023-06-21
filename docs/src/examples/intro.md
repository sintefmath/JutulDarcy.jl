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

In addition, there is experimental 3D visualization and well plotting found as a conditional extension. To get access to these features, you need at least Julia 1.9 and the `GLMakie` package loaded.
```julia
Pkg.add("GLMakie")
using GLMakie
```
