# High-level API

## Setup

### Meshes

JutulDarcy can use meshes that supported by Jutul. This includes the Cartesian and Unstructured meshes as well as any meshes in the more general [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) package.

```@docs
Jutul.CartesianMesh
Jutul.UnstructuredMesh
```

### Reservoir

```@docs
reservoir_domain
```

### Wells

```@docs
setup_well
setup_vertical_well
```

### Model

```@docs
setup_reservoir_model
```

### Initial state

```@docs
setup_reservoir_state
```

TODO: Write about hydrostatic equilbriation.

## Simulation

```@docs
simulate_reservoir
setup_reservoir_simulator
ReservoirSimResult
```
