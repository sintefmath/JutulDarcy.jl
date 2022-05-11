# JutulDarcy.jl
Darcy-scale and subsurface flow using [Jutul.jl](https://github.com/sintefmath/Jutul.jl) developed by the [Computational Geosciences group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/computational-geoscience/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Key features
- Written in pure Julia, with automatic differentiation
- High performance assembly and linear solvers, with support for two-stage CPR BILU(0)-CPR Krylov solvers
- Compositional, immiscible and black oil (experimental!) flow is supported and validated against existing simulators
- Unstructured grids and complex cases input from [the Matlab Reservoir Simulation Toolbox (MRST)](https://www.mrst.no).
- Support for general multisegment wells with rigorous mass balance, complex well limits and time-dependent controls
- 3D visualization of grids and wells
- Interactive plotting of well curves (somewhat experimental)

## A few of the packages used by Jutul and JutulDarcy
Jutul builds upon many of the excellent packages in the Julia ecosystem. Here are a few of them, and what they are used for:
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) implements the Dual number class used throughout the code
- [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) provides the iterative linear solvers
- [ILUZero.jl](https://github.com/mcovalt/ILUZero.jl/blob/master/src/ILUZero.jl) for ILU(0) preconditioners
- [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) for AMG preconditioners
- [Tullio.jl](https://github.com/mcabbott/Tullio.jl) for automatically optimized loops and [Polyester.jl]() for lightweight threads
- [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl) and [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) gives nice output to terminal 
- [Makie.jl](https://makie.juliaplots.org/) is used for the visualization features
- [MultiComponentFlash.jl](https://github.com/moyner/MultiComponentFlash.jl) provides many of the compositional features
- 
...and many more, both directly in the Project.toml file and indirectly!

## Getting started
As per Julia tradition, the documentation is mostly missing at the moment. This can be considered an early beta release. If you want to test the package, please start from one of the [examples found in this repository](https://github.com/sintefmath/JutulExamples.jl).

*Internals and undocumented functions are subject to change at this time. However, the main interface for the reservoir simulator itself seen in the examples should be fairly stable.*
