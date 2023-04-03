# Systems in JutulDarcy
JutulDarcy supports a number of different systems. These are [`JutulSystem`](@ref) instances that describe a particular type of physics for porous media flow. We describe these in roughly the order of complexity that they can model.

The general form of the flow systems we will discuss is a conservation law for ``N`` components on residual form:

```math
R = \frac{\partial}{\partial t} M_i + \nabla \cdot \vec{V}_i - Q_i, \quad \forall i \in \{1, \dots, N\}
```

Here, ``M_i`` is the conserved quantity (usually masses) for component ``i`` and ``\vec{V}_i`` the velocity of the conserved quantity. ``Q_i`` represents source terms that come from direct sources [`SourceTerm`](@ref), boundary conditions ([`FlowBoundaryCondition`](@ref)) or from wells ([`MultiSegmentWell`](@ref), [`SimpleWell`](@ref)).

### Implementation details
In the above the discrete version of ``M_i`` is implemented in the update function for [`TotalMasses`](@ref) that should by convention be named [`update_total_masses!`](@ref). The discrete component fluxes are implemented by [`component_mass_fluxes!`](@ref). The source terms are implemented by [`apply_forces_to_equation!`](@ref) for boundary conditions and sources, and [`update_cross_term_in_entity!`](@ref) for wells. We use Julia's multiple dispatch to pair the right implementation with the right physics system.

[``](@ref)
## Single-phase flow
The simplest form of porous media flow is the single-phase system.

```math
R(p) = \frac{\partial}{\partial t}( \rho \phi) + \nabla \cdot (\rho \vec{v}) - \rho q
```

``\rho`` is the phase mass density and ``\phi`` the apparent porosity of the medium, i.e. the void space in the rock available to flow. Where the velocity ``\vec{v}`` is given by Darcy's law that relates the the pressure gradient ``\nabla p`` and hydrostatic head to the velocity field:

```math
\vec{v} = - \frac{\mathbf{K}}{\mu} (\nabla p + \rho g \nabla z)
```

Here, ``\mathbf{K}`` is a positive-definite permeability tensor, ``\mu`` the fluid viscosity, ``g`` the magnitude of gravity oriented down and ``z`` the depth. 

!!! note "Single-phase implementation"
    The [`SinglePhaseSystem`](@ref) is a dedicated single phase system. This is mathematically equivalent to an [`ImmiscibleSystem`](@ref) when set up with a single phase. For single phase flow, the fluid [`Pressure`](@ref) is the primary variable in each cell. The equation supports two types of compressibility: That of the fluid where density is a function ``\rho(p)`` of pressure and that of the pores where the porosity ``\phi(p)`` changes with pressure.


!!! tip
    JutulDarcy uses the notion of depth rather than coordinate when defining buoyancy forces. This is consistent with the convention in the literature on subsurface flow.


## Multi-phase, immiscible flow
The flow systems immediately become more interesting if we add more phases. We can extend the above single-phase system by introducing the phase saturation of phase with label ``\alpha`` as ``S_\alpha``. The phase saturation represents the volumetric fraction of the rock void space occupied by the phase. If we consider a pair of phases ``\{n, w\}`` non-wetting and wetting we can write the system as

```math
R_\alpha = \frac{\partial}{\partial t} (S_\alpha \rho_\alpha \phi) + \nabla \cdot (\rho_\alpha \vec{v}_\alpha) - \rho_\alpha q_\alpha = 0, \quad \alpha \in \{n, w\}
```
This requires an additional closure such that the amount of saturation of all phases exactly fills the available fluid volume:
```math
S_w + S_n = 1, \quad 1 \ge S_\alpha \ge 0 \quad \alpha \in \{n, w\}
```
This equation is local and linear in the saturations and can be eliminated to produce the classical two-equation system for two-phase flow,
```math
R_n = \frac{\partial}{\partial t} ((1 - S_w) \rho_n \phi) + \nabla \cdot (\rho_n \vec{v}_n) - \rho_n q_n = 0,\\
R_w = \frac{\partial}{\partial t} (S_w \rho_w \phi) + \nabla \cdot (\rho_w \vec{v}_w) - \rho_w q_w = 0.
```
To complete this description we also need expressions for the phase fluxes. We use the standard multiphase extension of Darcy's law,
```math
\vec{v}_\alpha = - \mathbf{K} \frac{k_{r\alpha}}{\mu_\alpha} (\nabla p\alpha + \rho_\alpha g \nabla z)
```
Here, we have introduced the relative permeability of the phase ``k_{r\alpha}``, an empirical relationship between the saturation and the flow rate. Relative permeability is a complex topic with many different relationships and functional forms, but we limit the discussion to monotone, non-negative functions of their respective saturations, for example a simple Brooks-Corey type of ``k_{r\alpha}(S_\alpha) = S_\alpha^2``. We have also introduced separate phase pressures ``p_\alpha`` that account for capillary pressure, e.g. ``p_w = p_n + p_c(S_w)``.

!!! note "Immiscible implementation"
    The [`ImmiscibleSystem`](@ref) implements this system for any number of phases. The primary variables for this system is a single reference [`Pressure`](@ref) and phase [`Saturations`](@ref). As we do not solve for the volume closure equation, there is one less degree of freedom associated with the saturations than there are number of phases.
### ImmiscibleSystem

#### Primary variables
The 

## Multi-phase, pseudo-compositional flow

```math
R(p) = \frac{\partial}{\partial t}( \rho \phi) + \nabla \cdot (\rho \vec{v}) - q
```

## Multi-phase, multi-component equation-of-state flow
```math
R(p) = \frac{\partial}{\partial t}( \rho \phi) + \nabla \cdot (\rho \vec{v}) - q
```
### StandardBlackOilSystem
The [``](@ref)

#### Primary variables

### CompositionalSystem
The [``](@ref)


#### Primary variables

## Summary

## Thermal flow
Currently experimental.

# Solving the system

## Newton's method

## Linear solvers

### Single model (only porous medium)

### Multi model (porous medium with wells)
