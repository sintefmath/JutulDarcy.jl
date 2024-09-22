# JutulDarcy.jl - reservoir simulation and porous media flow in Julia

```@meta
CurrentModule = JutulDarcy
```

```@meta
DocTestSetup = quote
    using Jutul;
    using JutulDarcy;
end
```

Documentation for [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl). The documentation consists of a mixture of technical documentation, docstrings and examples. The examples are sorted by complexity. We suggest you start with [Gravity segregation example](Gravity segregation example).

To generate a folder that contains the examples locally, you can run the following code to create a folder `jutuldarcy_examples` in your current working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples()
```

Alternatively, a folder can be specified if you want the examples to be placed outside your present working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples("/home/username/")
```

JutulDarcy builds upon the general features found in [Jutul.jl](https://github.com/sintefmath/Jutul.jl). You may also find it useful to look at the [Jutul.jl documentation](https://sintefmath.github.io/Jutul.jl/dev/).

## Package docstring

```@docs
JutulDarcy
```

## Index

```@index
```
