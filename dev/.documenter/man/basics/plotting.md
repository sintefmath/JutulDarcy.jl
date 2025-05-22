
# Plotting and visualization {#Plotting-and-visualization}
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.plot_reservoir' href='#JutulDarcy.plot_reservoir'><span class="jlbinding">JutulDarcy.plot_reservoir</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_reservoir(model, states=missing; well_fontsize = 18, well_linewidth = 3, kwarg...)
```


Launch interactive plotter of reservoir + well trajectories in reservoir. Requires GLMakie.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/ext/ext_makie.jl#L64-L68" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.plot_well_results' href='#JutulDarcy.plot_well_results'><span class="jlbinding">JutulDarcy.plot_well_results</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_well_results(wr::WellResults)
plot_well_results(v::Vector{WellResults})
```


Launch interactive viewer for well results. Needs GLMakie to be loaded.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/ext/ext_makie.jl#L10-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.plot_reservoir_measurables' href='#JutulDarcy.plot_reservoir_measurables'><span class="jlbinding">JutulDarcy.plot_reservoir_measurables</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_reservoir_measurables(case::JutulCase, result::ReservoirSimResult)
```


Launch interactive viewer for reservoir measurables. Needs GLMakie to be loaded.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/ext/ext_makie.jl#L22-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.plot_reservoir_simulation_result' href='#JutulDarcy.plot_reservoir_simulation_result'><span class="jlbinding">JutulDarcy.plot_reservoir_simulation_result</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult; wells = true, reservoir = true)
```


Plot a reservoir simulation result. If `wells=true` well curves will be shown interactively. If `reservoir=true` the reservoir quantities will be visualized in 3D. These options can be combined.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/ext/ext_makie.jl#L32-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.plot_well!' href='#JutulDarcy.plot_well!'><span class="jlbinding">JutulDarcy.plot_well!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_well!(ax, mesh, w; color = :darkred)
```


Plot a given well that exists in mesh in Axis.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/ext/ext_makie.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

