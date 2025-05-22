
# Primary variables {#Primary-variables}

## Fluid systems {#Fluid-systems}

### General {#General}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.Pressure' href='#JutulDarcy.Pressure'><span class="jlbinding">JutulDarcy.Pressure</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Pressure(; max_abs = nothing, max_rel = 0.2, scale = si_unit(:bar), maximum = Inf, minimum = DEFAULT_MINIMUM_PRESSURE)
```


Pressure variable definition. `max_abs`/`max_rel` maximum allowable absolute/relative change over a Newton iteration, `scale` is a &quot;typical&quot; value used to regularize the linear system, `maximum` the largest possible value and `minimum` the smallest.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multiphase.jl#L83-L90" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ImmiscibleSaturation' href='#JutulDarcy.ImmiscibleSaturation'><span class="jlbinding">JutulDarcy.ImmiscibleSaturation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



A single saturation variable that represents the &quot;other&quot; phase in a three phase compositional system where two phases are predicted by an EoS


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multicomponent/variables/primary.jl#L69-L72" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Immiscible flow {#Immiscible-flow}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.Saturations' href='#JutulDarcy.Saturations'><span class="jlbinding">JutulDarcy.Saturations</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Saturations(; ds_max = 0.2)
```


Saturations as primary variable. `ds_max` controls maximum allowable saturation change between two successive Newton iterations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multiphase.jl#L105-L110" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Black-oil flow {#Black-oil-flow}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.BlackOilUnknown' href='#JutulDarcy.BlackOilUnknown'><span class="jlbinding">JutulDarcy.BlackOilUnknown</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BlackOilUnknown(dr_max = Inf, ds_max = Inf)
```


Variable defining the variable switching black-oil variable. The `ds_max` argument limits the maximum saturation change over a single Newton iteration when both a oileic and gaseous phase is present. The `dr_max` limits the maximum change in the undersaturated variable, taken relative to the maximum value of the undersaturated variable.

During simulation, this variable can take on the following interpretations: Gas saturation, Rs or Rv, depending on the phase conditions and miscibility model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/blackoil/variables/variables.jl#L71-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.BlackOilX' href='#JutulDarcy.BlackOilX'><span class="jlbinding">JutulDarcy.BlackOilX</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BlackOilX(sys::BlackOilVariableSwitchingSystem, p; sw = 0.0, so = 0.0, sg = 0.0, rs = 0.0, rv = 0.0, region = 1)
```


High level initializer for the black oil unknown degree of freedom. Will try to fill in the gaps unless system is really underspecified.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/blackoil/variables/variables.jl#L120-L125" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Compositional flow {#Compositional-flow}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.OverallMoleFractions' href='#JutulDarcy.OverallMoleFractions'><span class="jlbinding">JutulDarcy.OverallMoleFractions</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
OverallMoleFractions(;dz_max = 0.2)
```


Overall mole fractions definition for compositional. `dz_max` is the maximum allowable change in any composition during a single Newton iteration.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/multicomponent/variables/primary.jl#L23-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Wells {#Wells}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TotalMassFlux' href='#JutulDarcy.TotalMassFlux'><span class="jlbinding">JutulDarcy.TotalMassFlux</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TotalMassFlux(scale = si_unit(:day), max_abs = nothing, max_rel = nothing)
```


Variable normally used as primary variable. Represents the total mass flux going through a face. The typical usage is the mass flow through a segment of a [`MultiSegmentWell`](/man/basics/wells#JutulDarcy.MultiSegmentWell).

Note that the flow direction can often switch signs over a segment during a complex simulation. Setting `max_rel` to something other than `nothing` can therefore lead to severe convergence issues in the case of flow reversal.

**Fields (as keyword arguments)**
- `scale`: Scaling for variable
  
- `max_abs`: Max absolute change between Newton iterations
  
- `max_rel`: Maximum relative change between Newton iterations
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/wells.jl#L13-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## WellGroup / Facility {#WellGroup-/-Facility}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.TotalSurfaceMassRate' href='#JutulDarcy.TotalSurfaceMassRate'><span class="jlbinding">JutulDarcy.TotalSurfaceMassRate</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TotalSurfaceMassRate(max_absolute_change = nothing, max_relative_change = nothing)
```


Variable, typically representing the primary variable for a [`WellGroup`](/man/basics/wells#JutulDarcy.WellGroup). The variable is a single entry per well and solves for the total surface mass rate from a well to the facility model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/types.jl#L31-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

