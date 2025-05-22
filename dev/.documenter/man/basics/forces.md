
# Driving forces {#Driving-forces}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_reservoir_forces' href='#JutulDarcy.setup_reservoir_forces'><span class="jlbinding">JutulDarcy.setup_reservoir_forces</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_reservoir_forces(model; control = nothing, limits = nothing, set_default_limits = true, <keyword arguments>)
```


Set up driving forces for a reservoir model with wells


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L1344-L1348" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Source terms {#Source-terms}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.SourceTerm' href='#JutulDarcy.SourceTerm'><span class="jlbinding">JutulDarcy.SourceTerm</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SourceTerm(cell, value; fractional_flow = [1.0], type = MassSource)
```


Create source term in given `cell` with given total `value`.

The optional `fractional_flow` argument controls how this term is divided over components if used for inflow and should contain one entry per component in the system: (`number_of_components(system)`). `fractional_flow` should sum up to 1.0. The `type` argument should be an instance of the `FlowSourceType` enum, with interpretations as follows:
- `MassSource`: Source is directly interpreted as component masses.
  
- `StandardVolumeSource`: Source is volume at standard/surface conditions.  References densities are used to convert into mass sources.
  
- `VolumeSource`: Source is volume at in-situ / reservoir conditions.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/forces/sources.jl#L1-L17" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.FlowSourceType' href='#JutulDarcy.FlowSourceType'><span class="jlbinding">JutulDarcy.FlowSourceType</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



MassSource: Source is directly interpreted as component masses. StandardVolumeSource: Source is volume at standard/surface conditions. References densities are used to convert into mass sources. VolumeSource: Source is volume at in-situ / reservoir conditions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/types.jl#L340-L344" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Boundary conditions {#Boundary-conditions}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.FlowBoundaryCondition' href='#JutulDarcy.FlowBoundaryCondition'><span class="jlbinding">JutulDarcy.FlowBoundaryCondition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
FlowBoundaryCondition(
cell,
pressure = DEFAULT_MINIMUM_PRESSURE,
temperature = 298.15;
fractional_flow = nothing,
density = nothing,
trans_flow = 1e-12,
trans_thermal = 1e-6
)
```


Dirchlet boundary condition for constant values (pressure/temperature) at some inflow boundary


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/forces/bc.jl#L1-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.flow_boundary_condition' href='#JutulDarcy.flow_boundary_condition'><span class="jlbinding">JutulDarcy.flow_boundary_condition</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
flow_boundary_condition(cells, domain, pressures, temperatures = 298.15; kwarg...)
```


Add flow boundary conditions to a vector of `cells` for a given `domain` coming from `reservoir_domain`. The input arguments `pressures` and `temperatures` can either be scalars or one value per cell. Other keyword arguments are passed onto the `FlowBoundaryCondition` constructor.

The output of this function is a `Vector` of boundary conditions that can be passed on the form `forces = setup_reservoir_forces(model, bc = bc)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/forces/bc.jl#L84-L94" target="_blank" rel="noreferrer">source</a></Badge>

</details>

