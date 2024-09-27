# High-level API

## Setup

The basic outline of building a reservoir simulation problem consists of:

1. Making a mesh
2. Converting the mesh into a reservoir, adding properties
3. Add any number of wells
4. Setup a physical system and setup a reservoir model
5. Set up timesteps, well controls and other forces
6. Simulate!

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
