
# Wells and controls {#Wells-and-controls}

## Well setup routines {#Well-setup-routines}

Wells can be set up using the convenience functions [`setup_well`](/man/highlevel#JutulDarcy.setup_well) and [`setup_vertical_well`](/man/highlevel#JutulDarcy.setup_vertical_well). These routines act on the output from [`reservoir_domain`](/man/highlevel#JutulDarcy.reservoir_domain) and can set up both types of wells. We recommend that you use these functions instead of manually calling the well constructors.
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.WellGroup' href='#JutulDarcy.WellGroup'><span class="jlbinding">JutulDarcy.WellGroup</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
WellGroup(wells::Vector{Symbol}; can_shut_wells = true)
```


Create a well group that can control the given set of wells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L18-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Types of wells {#Types-of-wells}

### Simple wells {#Simple-wells}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SimpleWell' href='#JutulDarcy.SimpleWell'><span class="jlbinding">JutulDarcy.SimpleWell</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SimpleWell(reservoir_cells; <keyword arguments>)
```


Set up a simple well.

**Note**

[`setup_vertical_well`](/man/highlevel#JutulDarcy.setup_vertical_well) or [`setup_well`](/man/highlevel#JutulDarcy.setup_well) are the recommended way of setting up wells.

**Fields**
- `volume`
  
- `perforations`
  
- `surface`
  
- `name`
  
- `explicit_dp`
  
- `reference_depth`
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/types.jl#L400-L414" target="_blank" rel="noreferrer">source</a></Badge>

</details>


#### Equations {#Equations}

### Multisegment wells {#Multisegment-wells}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.WellDomain' href='#JutulDarcy.WellDomain'><span class="jlbinding">JutulDarcy.WellDomain</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract supertype for all well domains.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/types.jl#L371-L373" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.MultiSegmentWell' href='#JutulDarcy.MultiSegmentWell'><span class="jlbinding">JutulDarcy.MultiSegmentWell</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MultiSegmentWell(reservoir_cells, volumes, centers;
                N = nothing,
                name = :Well,
                perforation_cells = nothing,
                segment_models = nothing,
                segment_length = nothing,
                reference_depth = 0,
                dz = nothing,
                surface_conditions = default_surface_cond(),
                accumulator_volume = mean(volumes),
                )
```


Create well perforated in a vector of `reservoir_cells` with corresponding `volumes` and cell `centers`.

**Note**

[`setup_vertical_well`](/man/highlevel#JutulDarcy.setup_vertical_well) or [`setup_well`](/man/highlevel#JutulDarcy.setup_well) are the recommended way of setting up wells.

**Fields**
- `type`
  
- `volumes`
  
- `perforations`
  
- `neighborship`
  
- `top`
  
- `centers`
  
- `surface`
  
- `name`
  
- `segment_models`
  
- `material_thermal_conductivity`
  
- `material_density`
  
- `material_heat_capacity`
  
- `void_fraction`
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/types.jl#L460-L485" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SegmentWellBoreFrictionHB' href='#JutulDarcy.SegmentWellBoreFrictionHB'><span class="jlbinding">JutulDarcy.SegmentWellBoreFrictionHB</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Hagedorn and Brown well bore friction model for a segment.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/mswells_equations.jl#L2-L4" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.PotentialDropBalanceWell' href='#JutulDarcy.PotentialDropBalanceWell'><span class="jlbinding">JutulDarcy.PotentialDropBalanceWell</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PotentialDropBalanceWell(discretization)
```


Equation for the pressure drop equation in a multisegment well. This equation lives on the segment between each node and balances the pressure difference across the segment with the hydrostatic difference and well friction present in the current flow regime.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/mswells_equations.jl#L67-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.MixedWellSegmentFlow' href='#JutulDarcy.MixedWellSegmentFlow'><span class="jlbinding">JutulDarcy.MixedWellSegmentFlow</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Two point approximation with flux for wells


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/wells.jl#L5-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Well controls and limits {#Well-controls-and-limits}

### Types of well controls {#Types-of-well-controls}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.InjectorControl' href='#JutulDarcy.InjectorControl'><span class="jlbinding">JutulDarcy.InjectorControl</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
InjectorControl(target, mix, density = 1.0, phases = ((1, 1.0)), temperature = 293.15)
```


Well control that specifies injection into the reservoir. `target` specifies the type of target and `mix` defines the injection mass fractions for all species in the model during injection. 

For example, for a three-component system made up of CO2, H2O and H2, setting [0.1, 0.6, 0.3] would mean that the injection stream would contain 1 part CO2, 6 parts H2O and 3 parts H2 by mass. For an immiscible system (e.g. `LiquidPhase(), VaporPhase()`) the species corresponds to phases and [0.3, 0.7] would mean a 3 to 7 mixture of liquid and vapor by mass.

The density of the injected fluid at surface conditions is given by `density` which is defaulted to 1.0 if not given.

See also [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl), [`DisabledControl`](/man/basics/wells#JutulDarcy.DisabledControl).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L353-L368" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ProducerControl' href='#JutulDarcy.ProducerControl'><span class="jlbinding">JutulDarcy.ProducerControl</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ProducerControl(target)
```


Well control for production out of the reservoir. `target` specifies the type of target (for example `BottomHolePressureTarget()`).

See also [`DisabledControl`](/man/basics/wells#JutulDarcy.DisabledControl), [`InjectorControl`](/man/basics/wells#JutulDarcy.InjectorControl).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L440-L446" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.DisabledControl' href='#JutulDarcy.DisabledControl'><span class="jlbinding">JutulDarcy.DisabledControl</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DisabledControl()
```


Control that disables a well. If a well is disabled, it is disconnected from the surface network and no flow occurs between the well and the top side. Mass transfer can still occur inside the well, and between the well and the reservoir unless perforations are also closed by a [`PerforationMask`](/man/basics/wells#JutulDarcy.PerforationMask).

See also [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl), [`InjectorControl`](/man/basics/wells#JutulDarcy.InjectorControl).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L321-L330" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.replace_target' href='#JutulDarcy.replace_target'><span class="jlbinding">JutulDarcy.replace_target</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
replace_target(ctrl, new_target)
```


Create new well control using `ctrl` as a template that operates under `new_target`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L339-L343" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.default_limits' href='#JutulDarcy.default_limits'><span class="jlbinding">JutulDarcy.default_limits</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
default_limits(ctrl)
```


Create reasonable default limits for well control `ctrl`, for example to avoid BHP injectors turning into producers.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L306-L311" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Types of well targets {#Types-of-well-targets}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.BottomHolePressureTarget' href='#JutulDarcy.BottomHolePressureTarget'><span class="jlbinding">JutulDarcy.BottomHolePressureTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BottomHolePressureTarget(q, phase)
```


Bottom-hole pressure (bhp) target with target pressure value `bhp`. A well operating under a bhp constraint will keep the well pressure at the bottom hole (typically the top of the perforations) fixed at this value unless doing so would violate other constraints, like the well switching from injection to production when declared as an injector.

**Examples**

```julia
julia> BottomHolePressureTarget(100e5)
BottomHolePressureTarget with value 100.0 [bar]
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L98-L113" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SinglePhaseRateTarget' href='#JutulDarcy.SinglePhaseRateTarget'><span class="jlbinding">JutulDarcy.SinglePhaseRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SinglePhaseRateTarget(q, phase)
```


Single-phase well target with value `q` specified for `phase`.

**Examples**

```julia
julia> SinglePhaseRateTarget(0.001, LiquidPhase())
SinglePhaseRateTarget of 0.001 [m^3/s] for LiquidPhase()
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L119-L130" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SurfaceLiquidRateTarget' href='#JutulDarcy.SurfaceLiquidRateTarget'><span class="jlbinding">JutulDarcy.SurfaceLiquidRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SurfaceLiquidRateTarget(q)
```


Well target of specified liquid rate at surface conditions with value `q`. Typically used for a [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl) as you have full control over the mixture composition during injection.

Liquid rate, sometimes abbreviated LRAT, is made up of the phases that remain liquid at surface conditions. Typically, this will be water and oil if present in the model, but never different types of gas. If a producing becomes nearly or completely flooded by gas the well can go to very high or even infinite flows. It is therefore important to combine this control with a limit such as a bottom-hole-pressure constraint.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L137-L150" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SurfaceOilRateTarget' href='#JutulDarcy.SurfaceOilRateTarget'><span class="jlbinding">JutulDarcy.SurfaceOilRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SurfaceOilRateTarget(q)
```


Well target of specified oil rate with value `q` at surface conditions. Typically used for a [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl) as oil, for economic reasons, is rarely injected into the subsurface. Abbreviated as ORAT in some settings.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L163-L169" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SurfaceGasRateTarget' href='#JutulDarcy.SurfaceGasRateTarget'><span class="jlbinding">JutulDarcy.SurfaceGasRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SurfaceGasRateTarget(q)
```


Well target of specified gas rate with value `q` at surface conditions.

Often used for both [`InjectorControl`](/man/basics/wells#JutulDarcy.InjectorControl) [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl). Abbreviated as GRAT in some settings. If used for production it is important to also impose limits, as the well rate may become very high if there is little gas present.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L182-L191" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SurfaceWaterRateTarget' href='#JutulDarcy.SurfaceWaterRateTarget'><span class="jlbinding">JutulDarcy.SurfaceWaterRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SurfaceWaterRateTarget(q)
```


Well target of specified water rate with value `q` at surface conditions.

Often used for both [`InjectorControl`](/man/basics/wells#JutulDarcy.InjectorControl) [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl). If used for production it is important to also impose limits, as the well rate may become very high if there is little water present.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L204-L212" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TotalRateTarget' href='#JutulDarcy.TotalRateTarget'><span class="jlbinding">JutulDarcy.TotalRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TotalRateTarget(q)
```


Well target of specified total rate (sum of all phases) with value `q` at surface conditions.

Often used for both [`InjectorControl`](/man/basics/wells#JutulDarcy.InjectorControl) [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L225-L232" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.HistoricalReservoirVoidageTarget' href='#JutulDarcy.HistoricalReservoirVoidageTarget'><span class="jlbinding">JutulDarcy.HistoricalReservoirVoidageTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
HistoricalReservoirVoidageTarget(q, weights)
```


Historical RESV target for history matching cases. See [`ReservoirVoidageTarget`](/man/basics/wells#JutulDarcy.ReservoirVoidageTarget). For historical rates, the weights described in that target are computed based on the reservoir pressure and conditions at the previous time-step.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L263-L270" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ReservoirVoidageTarget' href='#JutulDarcy.ReservoirVoidageTarget'><span class="jlbinding">JutulDarcy.ReservoirVoidageTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ReservoirVoidageTarget(q, weights)
```


RESV target for history matching cases. The `weights` input should have one entry per phase (or pseudocomponent) in the system. The well control equation is then:

$|q_{ctrl} - \sum_i w_i q_i^s|$

where $q_i^s$ is the surface rate of phase $i$ and $w_i$ the weight of component stream $i$.

This constraint is typically set up from .DATA files for black-oil and immiscible cases.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L277-L290" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.DisabledTarget' href='#JutulDarcy.DisabledTarget'><span class="jlbinding">JutulDarcy.DisabledTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DisabledTarget(q)
```


Disabled target used when a well is under `DisabledControl()` only. The well will be disconnected from the surface.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L296-L301" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TotalReservoirRateTarget' href='#JutulDarcy.TotalReservoirRateTarget'><span class="jlbinding">JutulDarcy.TotalReservoirRateTarget</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TotalReservoirRateTarget(q)
```


Well target of specified total rate (sum of all phases) with value `q` at reservoir conditions.

Often used for both [`InjectorControl`](/man/basics/wells#JutulDarcy.InjectorControl) [`ProducerControl`](/man/basics/wells#JutulDarcy.ProducerControl).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L244-L251" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Implementation of well controls {#Implementation-of-well-controls}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.well_target' href='#JutulDarcy.well_target'><span class="jlbinding">JutulDarcy.well_target</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Well target contribution from well itself (disabled, zero value)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L266-L268" target="_blank" rel="noreferrer">source</a></Badge>



Well target contribution from well itself (bhp)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L273-L275" target="_blank" rel="noreferrer">source</a></Badge>



Well target contribution from well itself (surface volume, injector)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L280-L282" target="_blank" rel="noreferrer">source</a></Badge>



Well target contribution from well itself (reservoir volume, injector)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L296-L298" target="_blank" rel="noreferrer">source</a></Badge>



Well target contribution from well itself (surface volume, injector)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L309-L311" target="_blank" rel="noreferrer">source</a></Badge>



Well target contribution from well itself (surface volume, producer)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L322-L324" target="_blank" rel="noreferrer">source</a></Badge>



Well target contribution from well itself (RESV, producer)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/controls.jl#L358-L360" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Well outputs {#Well-outputs}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.print_well_result_table' href='#JutulDarcy.print_well_result_table'><span class="jlbinding">JutulDarcy.print_well_result_table</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
print_well_result_table(wr::WellResults, wells)
print_well_result_table(wr::WellResults, wells, outputs)
```


Print summary tables that show the well responses.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/well_results.jl#L152-L157" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Imposing limits on wells (multiple constraints) {#Imposing-limits-on-wells-multiple-constraints}

## Well forces {#Well-forces}

### Perforations and WI adjustments {#Perforations-and-WI-adjustments}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.PerforationMask' href='#JutulDarcy.PerforationMask'><span class="jlbinding">JutulDarcy.PerforationMask</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mask = PerforationMask(mask::Vector)
```


Create a perforation mask. This can be passed to [`setup_forces`](/ref/jutul#Jutul.setup_forces-Tuple{JutulModel}) for a well under the `mask` argument. The mask should equal the number of perforations in the well and is applied to the reference well indices in a multiplicative fashion. For example, if a well named `:Injector` has two perforations, the following mask would disable the first perforation and decrease the connection strength for the second perforation by 50%:

```julia
mask = PerforationMask([0.0, 0.5])
iforces = setup_forces(W, mask = mask)
forces = setup_reservoir_forces(model, control = controls, Injector = iforces)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L570-L584" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.Perforations' href='#JutulDarcy.Perforations'><span class="jlbinding">JutulDarcy.Perforations</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Perforations()
```


Entity that defines perforations: Connections from well cells to reservoir cells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L56-L61" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.compute_peaceman_index' href='#JutulDarcy.compute_peaceman_index'><span class="jlbinding">JutulDarcy.compute_peaceman_index</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_peaceman_index(g::T, K, r, pos; kwarg...) where T<:Jutul.JutulMesh
```


Compute the Peaceman index for a given mesh.

**Arguments**
- `g::JutulMesh`: Reservoir mesh
  
- `K`: Permeability tensor or scalar.
  
- `r`: Well radius.
  
- `pos`: Position of the well (index of cell or IJK truplet).
  
- `kwarg...`: Additional keyword arguments passed onto inner version of function.
  

**Returns**
- The computed Peaceman index.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/porousmedia.jl#L3-L19" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_peaceman_index(Δ, K, radius; kwargs...)
```


Compute the Peaceman well index for a given grid block.

**Arguments**
- `Δ`: The grid block size as a tuple `(dx, dy, dz)`
  
- `K`: The permeability of the medium (Matrix for full tensor, or scalar).
  
- `radius`: The well radius.
  

**Keyword Arguments**
- `dir::Symbol = :z`: Direction of the well, can be `:x`, `:y`, or `:z`.
  
- `net_to_gross = 1.0`: Net-to-gross ratio, used to scale the value for vertical directions.
  
- `constant = 0.14`: Constant used in the calculation of the equivalent radius. TPFA specific.
  
- `Kh = nothing`: Horizontal permeability, if not provided, it will be computed.
  
- `drainage_radius = nothing`: Drainage radius, if not provided, it will be computed.
  
- `skin = 0`: Skin factor, used to account for near-wellbore effects.
  
- `check = true`: Flag to check for negative well index values.
  

**Returns**
- `Float64`: The computed Peaceman well index.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/porousmedia.jl#L68-L90" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Other forces {#Other-forces}

Can use [`SourceTerm`](/man/basics/forces#JutulDarcy.SourceTerm) or [`FlowBoundaryCondition`](/man/basics/forces#JutulDarcy.FlowBoundaryCondition)
