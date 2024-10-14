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

Once a mesh has been set up, we can turn it into a reservoir with static properties:

```@docs
reservoir_domain
```

### Wells

Wells are most easily created using utilities that act directly on a reservoir domain:

```@docs
setup_well
setup_vertical_well
```

### Model

A single, option-heavy function is used to set up the reservoir model and default parameters:

```@docs
setup_reservoir_model
```

### Initial state

The initial state can be set up by explicitly setting all primary variables. JutulDarcy also contains functionality for initial hydrostatic equilibriation of the state, but this is at the moment most easily set up using input files with the `EQUIL` keyword.

```@docs
setup_reservoir_state
JutulDarcy.equilibriate_state
```

## Simulation

Simulating is done by either setting up a reservoir simulator and then simulating, or by using the convenience function that automatically sets up a simulator for you.

There are a number of different options available to tweak the tolerances, timestepping and behavior of the simulation. It is advised to read through the documentation in detail before running very large simulations.

```@docs
simulate_reservoir
setup_reservoir_simulator
ReservoirSimResult
```
