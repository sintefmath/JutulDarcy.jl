# Supported physical systems

JutulDarcy supports a number of different systems. These are [`JutulSystem`](https://sintefmath.github.io/Jutul.jl/dev/usage/#Jutul.JutulSystem) instances that describe a particular type of physics for porous media flow. We describe these in roughly the order of complexity that they can model.

## Summary

The general form of the flow systems we will discuss is a conservation law for ``N`` components on residual form:

```math
R = \frac{\partial}{\partial t} M_i + \nabla \cdot \vec{V}_i - Q_i, \quad \forall i \in \{1, \dots, N\}
```

Here, ``M_i`` is the conserved quantity (usually masses) for component ``i`` and ``\vec{V}_i`` the velocity of the conserved quantity. ``Q_i`` represents source terms that come from direct sources [`SourceTerm`](@ref), boundary conditions ([`FlowBoundaryCondition`](@ref)) or from wells ([`setup_well`](@ref), [`setup_vertical_well`](@ref)).

The following table gives an overview of the available features that are described in more detail below:

| System | Number of phases | Number of components | ``M`` | ``V`` |
|---|---|---|---|---|
| [`SinglePhaseSystem`](@ref) | 1 | 1 | ``\rho \phi`` | ``\rho \vec{v}`` |
| [`ImmiscibleSystem`](@ref) | Any | (Any) | ``S_\alpha \rho_\alpha \phi`` | ``\rho_\alpha \vec{v}_\alpha`` |
| [`StandardBlackOilSystem`](@ref) | 2-3 | (2-3) | ``\rho_o^s(b_g S_g + R_s b_o S_o)`` | ``b_g \vec{v}_g + R_s b_g \vec{v}_o`` |
| [`MultiPhaseCompositionalSystemLV`](@ref) | 2-3 | Any | ``\rho_l X_i S_l + \rho_v Y_i S_v`` | ``\rho_l X_i \vec{v}_l + \rho_v Y_i \vec{v}_v`` |

### Phases

Phases are defined using specific types. Some constructors take a list of phases present in the model. Phases do not contain any data themselves and the distinction between different phases applies primarily for well controls.

```@docs
LiquidPhase
AqueousPhase
VaporPhase
```

The phases are used by subtypes of the abstract superclass for multiphase flow systems:

```@docs
JutulDarcy.MultiPhaseSystem
```

### Implementation details

In the above the discrete version of ``M_i`` is implemented in the update function for [`JutulDarcy.TotalMasses`](@ref) that should by convention be named [`JutulDarcy.update_total_masses!`](@ref). The discrete component fluxes are implemented by [`JutulDarcy.component_mass_fluxes!`](@ref).

```@docs
JutulDarcy.TotalMasses
JutulDarcy.update_total_masses!
JutulDarcy.component_mass_fluxes!
```

The source terms are implemented by [`Jutul.apply_forces_to_equation!`](@ref) for boundary conditions and sources, and [`Jutul.update_cross_term_in_entity!`](@ref) for wells. We use Julia's multiple dispatch to pair the right implementation with the right physics system.

```@docs
Jutul.apply_forces_to_equation!
Jutul.update_cross_term_in_entity!
```

## Single-phase flow

```@docs
SinglePhaseSystem
```

The simplest form of porous media flow is the single-phase system.

```math
r(p) = \frac{\partial}{\partial t}( \rho \phi) + \nabla \cdot (\rho \vec{v}) - \rho q
```

``\rho`` is the phase mass density and ``\phi`` the apparent porosity of the medium, i.e. the void space in the rock available to flow. Where the velocity ``\vec{v}`` is given by Darcy's law that relates to the pressure gradient ``\nabla p`` and hydrostatic head to the velocity field:

```math
\vec{v} = - \frac{\mathbf{K}}{\mu} (\nabla p + \rho g \nabla z)
```

Here, ``\mathbf{K}`` is a positive-definite permeability tensor, ``\mu`` the fluid viscosity, ``g`` the magnitude of gravity oriented down and ``z`` the depth.

!!! note "Single-phase implementation"
    The [`SinglePhaseSystem`](@ref) is a dedicated single phase system. This is mathematically equivalent to an [`ImmiscibleSystem`](@ref) when set up with a single phase. For single phase flow, the fluid [`Pressure`](@ref) is the primary variable in each cell. The equation supports two types of compressibility: That of the fluid where density is a function ``\rho(p)`` of pressure and that of the pores where the porosity ``\phi(p)`` changes with pressure.

!!! tip
    JutulDarcy uses the notion of depth rather than coordinate when defining buoyancy forces. This is consistent with the convention in the literature on subsurface flow.

## Multi-phase, immiscible flow

```@docs
ImmiscibleSystem
```

The flow systems immediately become more interesting if we add more phases. We can extend the above single-phase system by introducing the phase saturation of phase with label ``\alpha`` as ``S_\alpha``. The phase saturation represents the volumetric fraction of the rock void space occupied by the phase. If we consider a pair of phases ``\{n, w\}`` non-wetting and wetting we can write the system as

```math
r_\alpha = \frac{\partial}{\partial t} (S_\alpha \rho_\alpha \phi) + \nabla \cdot (\rho_\alpha \vec{v}_\alpha) - \rho_\alpha q_\alpha = 0, \quad \alpha \in \{n, w\}
```

This requires an additional closure such that the amount of saturation of all phases exactly fills the available fluid volume:

```math
S_w + S_n = 1, \quad 1 \ge S_\alpha \ge 0 \quad \alpha \in \{n, w\}
```

This equation is local and linear in the saturations and can be eliminated to produce the classical two-equation system for two-phase flow,

```math
r_n = \frac{\partial}{\partial t} ((1 - S_w) \rho_n \phi) + \nabla \cdot (\rho_n \vec{v}_n) - \rho_n q_n = 0,\\
r_w = \frac{\partial}{\partial t} (S_w \rho_w \phi) + \nabla \cdot (\rho_w \vec{v}_w) - \rho_w q_w = 0.
```

To complete this description we also need expressions for the phase fluxes. We use the standard multiphase extension of Darcy's law,

```math
\vec{v}_\alpha = - \mathbf{K} \frac{k_{r\alpha}}{\mu_\alpha} (\nabla p\alpha + \rho_\alpha g \nabla z)
```

Here, we have introduced the relative permeability of the phase ``k_{r\alpha}``, an empirical relationship between the saturation and the flow rate. Relative permeability is a complex topic with many different relationships and functional forms, but we limit the discussion to monotone, non-negative functions of their respective saturations, for example a simple Brooks-Corey type of ``k_{r\alpha}(S_\alpha) = S_\alpha^2``. We have also introduced separate phase pressures ``p_\alpha`` that account for capillary pressure, e.g. ``p_w = p_n + p_c(S_w)``.

### Primary variables

!!! note "Immiscible implementation"
    The [`ImmiscibleSystem`](@ref) implements this system for any number of phases. The primary variables for this system is a single reference [`Pressure`](@ref) and phase [`Saturations`](@ref). As we do not solve for the volume closure equation, there is one less degree of freedom associated with the saturations than there are number of phases.

## Black-oil: Multi-phase, pseudo-compositional flow

```@docs
StandardBlackOilSystem
```

The black-oil equations is an extension of the immiscible description to handle limited miscibility between the phases. Originally developed for certain types of oil and gas simulation, these equations are useful when the number of components is low and tabulated values for dissolution and vaporization are available.

The assumptions of the black-oil model is that the "oil" and "gas" pseudo-components have uniform composition throughout the domain. JutulDarcy supports two- and three-phase black oil flow. The difference between two and three phases amounts to an additional immiscible aqueous phase that is identical to that of the previous section. For that reason, we focus on the miscible pseudo-components, oil:

```math
r_o = \rho_o^s \left( \frac{\partial}{\partial t}( (b_o S_o + R_v b_g (1 - S_o)) \phi) + \nabla \cdot ( b_o \vec{v}_o + R_v b_o \vec{v}_g) - q_o^s \right ),
```

and gas,

```math
r_g = \rho_g^s \left( \frac{\partial}{\partial t}( (b_g S_g + R_s b_o S_o) \phi) + \nabla \cdot ( b_g \vec{v}_g + R_s b_g \vec{v}_o) - q_g^s \right )
```

The model uses the notion of surface (or reference densities) ``\rho_o^s, \rho_g^s`` to define the densities of the component at specific pressure and temperature conditions where it is assumed that all "gas" has moved to the vapor phase and the defined "oil" is only found in the liquid phase. Keeping this definition in mind, the above equations can be divided by the surface densities to produce a surface volume balance equation where we have defined `b_o` and `b_g` as the dimensionless reciprocal formation volume factors that relate a volume at reservoir conditions to surface volumes and `R_s` for the dissolved volume of gas in the oil phase when brought to surface conditions. `R_v` is the same definition, but for oil vaporized into the gas phase.

!!! note "Blackoil implementation"
    The [`StandardBlackOilSystem`](@ref) implements the black-oil equations. It is possible to run cases with and without water, with and without ``R_s`` and with and without ``R_v``. The primary variables for the most general case is the reference [`Pressure`](@ref), an [`ImmiscibleSaturation`](@ref) for the aqueous phase and the special [`BlackOilUnknown`](@ref) that will represent either ``S_o``, ``R_s`` or ``R_v`` on a cell-by-cell basis depending on what phases are present and saturated.

A full description of the black-oil equations is outside the scope of this documentation. Please see [mrst-book-i](@cite) for more details.

## Compositional: Multi-phase, multi-component flow

```@docs
MultiPhaseCompositionalSystemLV
```

The more general case of multi-component flow is often referred to as a compositional model. The typical version of this model describes the fluid as a system of ``N`` components where the phases present and fluid properties are determined by an equation-of-state. This can be highly accurate if the equation-of-state is tuned for the mixtures that are encountered, but comes at a significant computational cost as the equation-of-state must be evaluated many times.

JutulDarcy implements a standard compositional model that assumes local instantaneous equilibrium and that the components are present in up to two phases with an optional immiscible phase added. This is sometimes referred to as a "simple water" or "dead water" description. By default the solvers use [MultiComponentFlash.jl](https://github.com/moyner/MultiComponentFlash.jl) to solve thermodynamic equilibrium. This package implements the generalized cubic approach and defaults to Peng-Robinson.

Assume that we have two phases liquid and vapor referred to as ``l`` and ``v`` with the Darcy flux given as in the preceeding sections. We can then write the residual equation for each of the ``M`` components by the liquid and vapor mole fractions ``X_i, Y_i```of that component as:

```math
r_i = \frac{\partial}{\partial t} \left( (\rho_l X_i S_l + \rho_v Y_i S_v) \phi \right) + \nabla \cdot (\rho_l X_i \vec{v}_l + \rho_v Y_i \vec{v}_v) - Q_i, \quad M \in \{1, \dots, M\}
```

For additional details, please see Chapter 8 - Compositional Simulation with the AD-OO Framework in [mrst-book-ii](@cite).

!!! note "Compositional implementation"
    The [`MultiPhaseCompositionalSystemLV`](@ref) implements the compositional model. The primary variables for the most general case is the reference [`Pressure`](@ref), an [`ImmiscibleSaturation`](@ref) for the optional immiscible phase and ``M-1`` [`OverallMoleFractions`](@ref).

## Multi-phase thermal flow

Thermal effects are modelled as an additional residual equation that comes in addition to the flow equations.

```math
r_i = \frac{\partial}{\partial t}\Bigl(\rho_r U_r (1-\phi) + \sum_\alpha \rho_\alpha S_\alpha U_\alpha \phi \Bigr) + \nabla \cdot \left (\sum_\alpha ( H_\alpha \rho_\alpha v_\alpha - S_\alpha \lambda_\alpha \nabla T)-\lambda_r \nabla T \right) - Q_e
```

```@docs
JutulDarcy.add_thermal_to_model!
```
