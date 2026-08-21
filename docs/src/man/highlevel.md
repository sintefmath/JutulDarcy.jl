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

JutulDarcy can use meshes that are supported by Jutul. This includes the Cartesian ([`Jutul.CartesianMesh`](@ref)) and Unstructured meshes ([`Jutul.UnstructuredMesh`](@ref)), meshes from Gmsh ([`Jutul.mesh_from_gmsh`](@ref)) and meshes from [MRST](https://www.mrst.no) ([`Jutul.MRSTWrapMesh`](@ref)) and different types of reservoir meshes parsed by [GeoEnergyIO.jl](https://sintefmath.github.io/GeoEnergyIO.jl/stable/) (RESQML, GRDECL/corner-point, TOPS/DX/DY/DZ and so on).

For convenience purposes, the `reservoir_mesh` function can be used to construct reservoir meshes of different types, and serves as an exensible interface that can be used by other packages to add new meshes to JutulDarcy. The same function can be used to retrieve meshes from already constructed objects.

```@docs
reservoir_mesh
```

### Reservoir

Once a mesh has been set up, we can turn it into a reservoir with static properties:

```@docs
reservoir_domain
get_1d_reservoir
```

### Fractures

To model fracture flow, we can construct a fracture domain from e.g, an embedded mesh [`Jutul.EmbeddedMeshes.EmbeddedMesh`](@ref) made up of a subset of the faces of a 3D reservoir mesh.

```@docs
fracture_domain
```

### Wells

Wells are most easily created using utilities that act directly on a reservoir domain:

```@docs
setup_well
setup_vertical_well
```

#### Setting up well controls

```@docs
setup_injector_control
setup_producer_control
```

### Defining forces

```@docs
setup_reservoir_forces
```

#### Source terms

```@docs
SourceTerm
JutulDarcy.FlowSourceType
```

#### Boundary conditions

```@docs
FlowBoundaryCondition
flow_boundary_condition
```

### Model

A single, option-heavy function is used to set up the reservoir model and default parameters:

```@docs
setup_reservoir_model
```

#### Fractured reservoirs

Similar functionality can be used to set up a fractured reservoir model using `setup_fractured_reservoir_model`. The function takes a fracture domain (constructed using `fracture_domain`)  as its second argument, in a addition to the keyword arguments of `setup_reservoir_model`

```@docs
setup_fractured_reservoir_model
```

### Properties and functions

#### Relative permeability

```@docs
set_relative_permeability!
JutulDarcy.StoneIIMethod
JutulDarcy.SaturationWeightedOilRelperm
JutulDarcy.StoneIMethod
```

#### Capillary pressure

```@docs
set_capillary_pressure!
```

### Initial state

The initial state can be set up by explicitly setting all primary variables. JutulDarcy also contains functionality for initial hydrostatic equilibriation of the state, which is either done by setting up `EquilibriumRegion` instances that are passed to `setup_reservoir_state`, or by using an input file with the `EQUIL` keyword.

```@docs
setup_reservoir_state
EquilibriumRegion
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
