
```@meta
CurrentModule = JutulDarcy
```

```@meta
DocTestSetup = quote
    using Jutul;
    using JutulDarcy;
end
```

# JutulDarcy

```@docs
JutulDarcy
```

Documentation for [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl). The documentation is currently limited to docstrings and a series of examples. The examples are sorted by complexity. We suggest you start with [Gravity segregation example](Gravity segregation example).

!!! info "Note about units"
    JutulDarcy does currently not make us of conversion factors or explicit
    units can in principle use any consistent unit system. Some default scaling
    of variables assume that the magnitude pressures and velocities roughly
    match that of strict SI (e.g. Pascals and cubic meters per second). These
    scaling factors are primarily used when iterative linear solvers are used.

JutulDarcy builds upon the general features found in [Jutul.jl](https://github.com/sintefmath/Jutul.jl). You may also find it useful to look at the [Jutul.jl documentation](https://sintefmath.github.io/Jutul.jl/dev/).

## Reading input files

It is also possible to read cases that have been set up in MRST (see [`setup_case_from_mrst`](@ref) and [`simulate_mrst_case`](@ref)) or from .DATA files (see [`parse_data_file`](@ref) and [`simulate_data_file`](@ref))

```@index
```
