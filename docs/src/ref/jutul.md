# Documentation from Jutul.jl

`JutulDarcy.jl` builds upon `Jutul.jl`, which takes care of the heavy lifting in terms of meshes, discretizations and solvers. You can use `JutulDarcy.jl` without
knowing the inner workings of `Jutul.jl`, but if you want to dive under the hood the [Jutul.jl manual](https://sintefmath.github.io/Jutul.jl/dev/) and [Jutul.jl docstrings](https://sintefmath.github.io/Jutul.jl/dev/docstrings/) may be useful.

We also include a few choice docstrings here that are extensively used in the code:
```@docs
Jutul.setup_forces
Jutul.JutulCase
Jutul.DataDomain
Jutul.setup_state
Jutul.simulator_config
Jutul.VectorVariables
Jutul.simulate!
Jutul.simulate
Jutul.Simulator
Jutul.SimulatorModel
Jutul.MultiModel
```
