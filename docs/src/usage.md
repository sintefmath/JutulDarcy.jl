# Systems in JutulDarcy
JutulDarcy supports a number of different systems. These are [`JutulSystem`](@ref) instances that describe a particular type of physics for porous media flow. We describe these in roughly the order of complexity that they can model.

The general form of the flow systems we will discuss is a conservation law for ``N`` components:

```math
\frac{\partial}{\partial t} M_i + \nabla \cdot \vec{V}_i - Q_i = 0, \quad \forall i \in \{1, \dots, N\}
```

Here, ``M_i`` is the conserved quantity (usually masses) for component ``i`` and ``\vec{V}_i`` the velocity of the conserved quantity. ``Q_i`` represents source terms that come from direct sources [`SourceTerm`](@ref), boundary conditions ([`FlowBoundaryCondition`](@ref)) or from wells ([`MultiSegmentWell`](@ref), [`SimpleWell`](@ref)).

### Implementation details
In the above the discrete version of ``M_i`` is implemented in the update function for [`TotalMasses`](@ref) that should by convention be named [`update_total_masses!`](@ref). The discrete component fluxes are implemented by [`component_mass_fluxes!`](@ref). The source terms are implemented by [`apply_forces_to_equation!`](@ref) for boundary conditions and sources, and [`update_cross_term_in_entity!`](@ref) for wells. We use Julia's multiple dispatch to pair the right implementation with the right physics system.

[``](@ref)
## Single-phase flow
The simplest form of porous media flow is the single-phase system.

```math
\frac{\partial}{\partial t}( \rho \phi) + \nabla \cdot (\rho \vec{v}) - q = 0
```

``\rho`` is the phase mass density and ``\phi`` the apparent porosity of the medium, i.e. the void space in the rock available to flow. Where the velocity ``\vec{v}`` is given by Darcy's law that relates the the pressure gradient ``\nabla p`` and hydrostatic head to the velocity field:

```math
\vec{v} = - \frac{\mathbf{K}}{\mu} (\nabla p + \rho g \nabla z)
```

Here, ``\mathbf{K}`` is a positive-definite permeability tensor, ``\mu`` the fluid viscosity, ``g`` the magnitude of gravity oriented down and ``z`` the depth. 


!!! note
    JutulDarcy uses the notion of depth rather than coordinate when defining buoyancy forces. This is consistent with the convention in the literature on subsurface flow.

### SinglePhaseSystem

The [`SinglePhaseSystem`](@ref) is a dedicated single phase system. This is mathematically equivalent to an [`ImmiscibleSystem`](@ref) when set up with a single phase.
#### Primary variables
For single phase flow, the fluid [`Pressure`](@ref) is the primary variable in each cell. The equation supports two types of compressibility: That of the fluid where density is a function ``\rho(p)`` of pressure and that of the pores where the porosity ``\phi(p)`` changes with pressure.

## Multi-phase flow

### ImmiscibleSystem
The [``](@ref)


!!! tip 
    Immiscible systems can be created with a single phase and are then equivialent in behavior to the single phase system.


## Multi-phase, multi-component flow

### StandardBlackOilSystem
The [``](@ref)

### CompositionalSystem
The [``](@ref)

## Summary

## Thermal flow
Currently experimental.