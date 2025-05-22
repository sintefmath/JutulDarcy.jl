
# Parameters {#Parameters}

## General {#General}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.FluidVolume' href='#JutulDarcy.FluidVolume'><span class="jlbinding">JutulDarcy.FluidVolume</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
FluidVolume()
```


Variable typically taken to be a parameter. Represents the per-cell volume that where multiphase flow can occur. For a well, this is the volume inside the well-bore free flow can occur. For a porous medium, this is the void space inside the pores that is, to some extent, connected and open to flow (effective pore-volume).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/pvt.jl#L125-L133" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Reservoir parameters {#Reservoir-parameters}

### Transmissibility {#Transmissibility}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.Transmissibilities' href='#JutulDarcy.Transmissibilities'><span class="jlbinding">JutulDarcy.Transmissibilities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Transmissibilities()
```


Variable/parameter used to define the cell-to-cell transmissibilities when using a two-point flux approximation scheme.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multiphase.jl#L183-L188" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.reservoir_transmissibility' href='#JutulDarcy.reservoir_transmissibility'><span class="jlbinding">JutulDarcy.reservoir_transmissibility</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
reservoir_transmissibility(d::DataDomain)
reservoir_transmissibility(d::DataDomain; version = :xyz)
```


Special transmissibility function for reservoir simulation that handles additional complexity present in industry grade models such as fault multipliers, net-to-gross etc. The input should be a `DataDomain` instance returned from [`reservoir_domain`](/man/highlevel#JutulDarcy.reservoir_domain)

The keyword argument `version` can be `:xyz` for permeability tensors that are aligned with coordinate directions or `:ijk` to interpreted the permeability as a diagonal tensor aligned with the logical grid directions. The latter choice is only meaningful for a diagonal tensor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L1825-L1838" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Other {#Other}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TwoPointGravityDifference' href='#JutulDarcy.TwoPointGravityDifference'><span class="jlbinding">JutulDarcy.TwoPointGravityDifference</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TwoPointGravityDifference()
```


Parameter representing the difference in gravity on an instance of `Faces` between two `Cells`. If the phase flux is written as

$v = - K \nabla (p + \rho g \nabla z)$

this term represent the discretized analogue of $\rho g \nabla z$.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multiphase.jl#L269-L278" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ConnateWater' href='#JutulDarcy.ConnateWater'><span class="jlbinding">JutulDarcy.ConnateWater</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConnateWater()
```


Parameter for connate water per cell. Used in some three-phase relative permeability evaluators.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multiphase.jl#L132-L137" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.EndPointScalingCoefficients' href='#JutulDarcy.EndPointScalingCoefficients'><span class="jlbinding">JutulDarcy.EndPointScalingCoefficients</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Type that defines endpoint scaling parameters for relative permeabilities (either drainage or imbibiton).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/variables/relperm/endscale.jl#L78-L81" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Well parameters {#Well-parameters}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.WellIndices' href='#JutulDarcy.WellIndices'><span class="jlbinding">JutulDarcy.WellIndices</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
WellIndices()
```


Parameter for the connection strength between a well and the reservoir for a given perforation. Typical values come from a combination of Peaceman&#39;s formula, upscaling and/or history matching.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L64-L70" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.PerforationGravityDifference' href='#JutulDarcy.PerforationGravityDifference'><span class="jlbinding">JutulDarcy.PerforationGravityDifference</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PerforationGravityDifference()
```


Parameter for the height difference from the wellbore and the connected node in the well.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L82-L87" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Thermal {#Thermal}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.FluidThermalConductivities' href='#JutulDarcy.FluidThermalConductivities'><span class="jlbinding">JutulDarcy.FluidThermalConductivities</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
FluidThermalConductivities()
```


Variable defining the fluid component conductivity.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/thermal/thermal.jl#L179-L183" target="_blank" rel="noreferrer">source</a></Badge>

</details>

