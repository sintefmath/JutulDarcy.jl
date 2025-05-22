
# High-level API {#High-level-API}

## Setup {#Setup}

The basic outline of building a reservoir simulation problem consists of:
1. Making a mesh
  
2. Converting the mesh into a reservoir, adding properties
  
3. Add any number of wells
  
4. Setup a physical system and setup a reservoir model
  
5. Set up timesteps, well controls and other forces
  
6. Simulate!
  

### Meshes {#Meshes}

JutulDarcy can use meshes that are supported by Jutul. This includes the Cartesian ([`Jutul.CartesianMesh`](/ref/jutul#Jutul.CartesianMesh)) and Unstructured meshes ([`Jutul.CartesianMesh`](/ref/jutul#Jutul.CartesianMesh)), meshes from Gmsh ([`Jutul.mesh_from_gmsh`](/ref/jutul#Jutul.mesh_from_gmsh)), meshes from [MRST](https://www.mrst.no) ([`Jutul.MRSTWrapMesh`](/man/basics/input_files#Jutul.MRSTWrapMesh)), and meshes from the [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) package.

### Reservoir {#Reservoir}

Once a mesh has been set up, we can turn it into a reservoir with static properties:
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.reservoir_domain' href='#JutulDarcy.reservoir_domain'><span class="jlbinding">JutulDarcy.reservoir_domain</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
reservoir_domain(g; permeability = convert_to_si(0.1, :darcy), porosity = 0.1, kwarg...)
```


Set up a `DataDomain` instance for given mesh or other representation `g`. `permeability` and `porosity` are then added to the domain. If scalars are passed, they are expanded to cover all cells. Arrays are asserted to match all cells. Permeability is either one value per cell (diagonal scalar), one value per dimension given in each row (for a diagonal tensor) or a vector that represents a compact full tensor representation (6 elements in 3D, 3 in 2D).

**Default data and their values**

|                         Name |                                Explanation |       Unit | Default |
| ----------------------------:| ------------------------------------------:| ----------:| -------:|
|               `permeability` |         Rock ability to conduct fluid flow |      $m^2$ |  100 mD |
|                   `porosity` |   Rock void fraction open to flow (0 to 1) |          - |     0.3 |
|               `rock_density` |                       Mass density of rock | $kg^3/m^3$ |  2000.0 |
|         `rock_heat_capacity` |             Specific heat capacity of rock | $J/(kg K)$ |   900.0 |
|  `rock_thermal_conductivity` |                  Heat conductivity of rock |    $W/m K$ |     3.0 |
| `fluid_thermal_conductivity` |          Heat conductivity of fluid phases |    $W/m K$ |     0.6 |
|    `component_heat_capacity` | Specific heat capacity of fluid components | $J/(kg K)$ |  4184.0 |


Note that the default values are taken to be roughly those of water for fluid phases and sandstone for those of rock. Choice of values can severely impact your simulation results - take care to check the values that your physical system makes use of!

**Optional values**

|                          Name |                               Explanation |    Unit | Default |
| -----------------------------:| -----------------------------------------:| -------:| -------:|
|                `net_to_gross` |   Magnitude of porosity available to flow |       - |     1.0 |
|                   `diffusion` |  Diffusion coefficient for each component | $m^2/s$ |     0.0 |
|   `transmissibility_override` |   Transmissibility override for each face |   $m^2$ |     NaN |
| `transmissibility_multiplier` | Transmissibility multiplier for each face |       - |     1.0 |


These values are optional and will only be added if specified. For e.g. net-to-gross and transmissibility multipliers a default value of 1.0 will be assumed in the code if it is not present.

The transmissibility override can be used to override the transmissibility calculated from geometry and other properties. This is useful for example if you have an externally computed model. The values must be given for every face on the mesh. Values that are NaN or Inf will be treated as missing and the standard transmissibility calculator will be used instead for those faces.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L48-L92" target="_blank" rel="noreferrer">source</a></Badge>



```julia
reservoir_domain(m::Union{SimulationModel, MultiModel})
```


Get reservoir domain embedded in model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L153-L157" target="_blank" rel="noreferrer">source</a></Badge>



```julia
reservoir_domain(case::JutulCase)
```


Get reservoir domain from a reservoir simulation case.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L163-L167" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.get_1d_reservoir' href='#JutulDarcy.get_1d_reservoir'><span class="jlbinding">JutulDarcy.get_1d_reservoir</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_1d_reservoir(nc;
    L = 1.0,
    perm = 9.8692e-14, # 0.1 darcy
    poro = 0.1,
    area = 1.0,
    z_max = nothing
)
```


Utility function for setting up a 1D reservoir domain with `nc` cells and length `L`. The [`reservoir_domain`](/man/highlevel#JutulDarcy.reservoir_domain) function is generally preferred and this function is kept for backwards compatibility.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/porousmedia_grids.jl#L32-L44" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Wells {#Wells}

Wells are most easily created using utilities that act directly on a reservoir domain:
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_well' href='#JutulDarcy.setup_well'><span class="jlbinding">JutulDarcy.setup_well</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_well(D::DataDomain, reservoir_cells; skin = 0.0, Kh = nothing, radius = 0.1, dir = :z, name = :Well)
w = setup_well(D, 1, name = :MyWell)         # Cell 1 in the grid
w = setup_well(D, (2, 5, 1), name = :MyWell) # Position (2, 5, 1) in logically structured mesh
w2 = setup_well(D, [1, 2, 3], name = :MyOtherWell)
```


Set up a well in `reservoir_cells` with given skin factor and radius. The order of cells matter as it is treated as a trajectory.

The `name` keyword argument can be left defaulted if your model will only have a single well (named `:Well`). It is highly recommended to provide this whenever a well is set up.

`reservoir_cells` can be one of the following: A Vector of cells, a single cell, a Vector of `(I, J, K)` Tuples or a single Tuple of the same type.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/wells.jl#L71-L87" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_vertical_well' href='#JutulDarcy.setup_vertical_well'><span class="jlbinding">JutulDarcy.setup_vertical_well</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_vertical_well(D::DataDomain, i, j; name = :MyWell, <kwarg>)
```


Set up a vertical well with a [`DataDomain`](/ref/jutul#Jutul.DataDomain-Tuple{JutulDomain}) input that represents the porous medium / reservoir where the wells it to be placed. See [`SimpleWell`](/man/basics/wells#JutulDarcy.SimpleWell), [`MultiSegmentWell`](/man/basics/wells#JutulDarcy.MultiSegmentWell) and [`setup_well`](/man/highlevel#JutulDarcy.setup_well) for more details about possible keyword arguments.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/wells.jl#L273-L280" target="_blank" rel="noreferrer">source</a></Badge>



```julia
setup_vertical_well(g, K, i, j; heel = 1, toe = grid_dims_ijk(g)[3], kwarg...)
```


Set up a vertical well for given grid `g` and permeability `K` at logical indices `i, j` perforating all cells starting at k-logical index `heel` to `toe`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/facility/wells/wells.jl#L294-L300" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Model {#Model}

A single, option-heavy function is used to set up the reservoir model and default parameters:
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_reservoir_model' href='#JutulDarcy.setup_reservoir_model'><span class="jlbinding">JutulDarcy.setup_reservoir_model</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
model, parameters = setup_reservoir_model(reservoir, system; wells = [], <keyword arguments>)
model, parameters = setup_reservoir_model(reservoir, system; wells = [w1, w2], backend = :csr, <keyword arguments>)
```


Set up a reservoir `MultiModel` for a given reservoir `DataDomain` typically set up from  [`reservoir_domain`](/man/highlevel#JutulDarcy.reservoir_domain) and an optional vector of wells that are created using [`setup_vertical_well`](/man/highlevel#JutulDarcy.setup_vertical_well) and  [`setup_well`](/man/highlevel#JutulDarcy.setup_well).

The routine automatically sets up a facility and couples the wells with the reservoir and that facility.

**Keyword arguments**

**Basic model setup**
- `wells=[]`: Vector of wells (e.g. from [`setup_well`](/man/highlevel#JutulDarcy.setup_well)) that are to be used in the model. Each well must have a unique name.
  
- `extra_out=true`: Return both the model and the parameters instead of just the model.
  
- `thermal = false`: Add additional equations for conservation of energy and temperature as a primary variable.
  
- `kgrad=nothing`: Type of spatial discretization to use:
  - `:tpfa` or `nothing` gives standard two-point flux approximation (TPFA) with hard-coded two-point assembly
    
  - `:tpfa_test` gives TPFA with specialized finite-volume assembly. Should be similar in performance to `:tpfa`, but does not make use of threads.
    
  - `:avgmpfa` gives a consistent linear MPFA scheme that is more accurate for meshes with anisotropic perm or non-orthogonal cells than `:tpfa`.
    
  - `:ntpfa` gives a consistent nonlinear MPFA scheme (nonlinear version of `:avgmpfa` that preserves monotonicity)
    
  
- `upwind=nothing`: Type of upwinding to use. Can be `:spu` or `nothing` for standard upwinding or `:weno` for a second-order weighted essentially non-oscillatory scheme.
  
- `extra_outputs=Symbol[]`: Extra output variables for reservoir model. Defaults to &quot;typical&quot; values seen in reservoir simulation. Valid values: Vector of symbols to be output, `true` for all variables and `false` for the minimal set required to restart simulations (typically only the primary variables and mass of each component)
  

**Advanced model setup**

Advanced options govern internals of the simulator, like type of automatic differentation, how equations are linearized and so on. These should not impact simulation results beyond what is allowed for the model tolerances, but can impact simulation speed.
- `split_wells=false`: Add a facility model for each well instead of one facility model that controls all wells. This must be set to `true` if you want to use MPI or nonlinear domain decomposition.
  
- `backend=:csr`: Backend to use. Can be `:csc` for serial compressed sparse column CSC matrix, `:csr` for parallel compressed sparse row matrix. `:csr` is a bit faster and is recommended when using MPI, HYPRE or multiple threads. `:csc` uses the default Julia format and is interoperable with other Julia libraries.
  
- `context=DefaultContext()`: Context used for entire model. Not recommended to set up manually, use `backend` instead.
  
- `assemble_wells_together=true`: Assemble wells in a single big matrix rather than many small matrices.
  
- `block_backend=true`: Use block sparse representation. This is needed by the iterative solvers and corresponding preconditioners. Setting this to `false` will result in a direct solver being used. In addition, equations will be assembled in an order similar to that of MRST (equation major instead of cell major).
  
- `general_ad=false`: Use more general form of AD. Will result in slower execution speed than if set to true, but can be useful when working with custom discretizations.
  
- `discretization_arg=NamedTuple()`: Additional keyword arguments passed onto `discretized_domain_tpfv_flow` when setting up discretizations.
  

**Increment and variable options**

These options govern the range of values and the maximum allowable change of properties over a single Newton iteration. Changing values for maximum change will not change the equations themselves, but the values will change the rate of nonlinear solver convergence. Typically, smaller values are more conservative and reduce numerical difficulties, but can significantly increase the number of iterations and the reduce the length of the average time-step. Setting very small values can make it infeasible to solve the problems in a reasonable time.

Note that relative values are usually given relative to the cell value. If your expected output values are close to zero (e.g. for near-atmospheric pressure) low values can lead to slow convergence.
- `dp_max_abs=nothing`: Maximum allowable pressure change in SI units (Pa)
  
- `dp_max_rel=0.2`: Maximum allowable relative pressure change (default is 20%)
  
- `dp_max_abs_well=convert_to_si(50, :bar)`: Maximum allowable pressure change for wells in SI units (Pa)
  
- `dp_max_rel_well=nothing`: Maximum allowable relative pressure change in well
  
- `ds_max=0.2`: Maximum change in saturations
  
- `dz_max=0.2`: Maximum change in composition (for compositional models only)
  
- `p_min=JutulDarcy.DEFAULT_MINIMUM_PRESSURE`: Minimum pressure in model (hard limit)
  
- `p_max=Inf`: Maximum pressure in model (hard limit)
  
- `dr_max=Inf`: Maximum change in Rs/Rv for blackoil models over a Newton iteration. Taken relative to the saturated value of the cell.
  
- `dT_max_rel=nothing`: Maximum relative change in temperature (JutulDarcy uses Kelvin, so comments about changing limits near zero above does not apply to typical reservoir temperatures)
  
- `dT_max_abs=50.0`: Maximum absolute change in temperature (in °K/°C)
  
- `T_min=convert_to_si(0.0, :Celsius)`: Minimum temperature in model (hard limit)
  
- `fast_flash=false`: Shorthand to enable `flash_reuse_guess` and `flash_stability_bypass`. These options can together speed up the time spent in flash solver for compositional models. Options are based on &quot;Increasing the Computational Speed of Flash Calculations With Applications for Compositional, Transient Simulations&quot; by Rasmussen et al (2006).
  
- `flash_reuse_guess=fast_flash`: Reuse previous flash guess when a cell remains in two-phase.
  
- `flash_stability_bypass=fast_flash`: Bypass stability testing for cells outside the two-phase and shadow region.
  
- `can_shut_wells=true`: Configure facility to allow shutting wells that repeatedly get rates with the wrong side. Disabling this can make certain models infeasible to simulate, but it can be useful to do so for simple models where you know that the wells should be operational.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L266-L378" target="_blank" rel="noreferrer">source</a></Badge>



```julia
setup_reservoir_model(reservoir::DataDomain, model_template::MultiModel; wells = [])
setup_reservoir_model(model_template::MultiModel; wells = [])
```


Set up a reservoir model with another model as a template. The template model is used to define the parameters and variables so that the resulting model is as similar to the original model as possible. The main purpose of this model is to &quot;resetup&quot; a model with for example a new set of wells.

It is also possible to pass a `reservoir` domain to set up the model with a new domain and new wells, copying over properties and secondary variables.

Note that the transfer process is not perfect and some variables might not be copied correctly if you are using a highly customized model. For instance, the treatment of regions is quite simple as it is based on the field name `region` that is initialized to one for each cell.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L587-L603" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Initial state {#Initial-state}

The initial state can be set up by explicitly setting all primary variables. JutulDarcy also contains functionality for initial hydrostatic equilibriation of the state, which is either done by setting up `EquilibriumRegion` instances that are passed to `setup_reservoir_state`, or by using an input file with the `EQUIL` keyword.
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_reservoir_state' href='#JutulDarcy.setup_reservoir_state'><span class="jlbinding">JutulDarcy.setup_reservoir_state</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_reservoir_state(model, <keyword arguments>)
# Ex: For immiscible two-phase
setup_reservoir_state(model, Pressure = 1e5, Saturations = [0.2, 0.8])
```


Convenience constructor that initializes a state for a `MultiModel` set up using [`setup_reservoir_model`](/man/highlevel#JutulDarcy.setup_reservoir_model). The main convenience over [`setup_state`](/ref/jutul#Jutul.setup_state-Tuple{JutulModel,%20Vararg{Any}}) is only the reservoir initialization values need be provided: wells are automatically initialized from the connected reservoir cells.

As an alternative to passing keyword arguments, a `Dict{Symbol, Any}` instance can be sent in as a second, non-keyword argument.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/state0.jl#L167-L179" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.EquilibriumRegion' href='#JutulDarcy.EquilibriumRegion'><span class="jlbinding">JutulDarcy.EquilibriumRegion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
EquilibriumRegion(model, p_datum = missing,
    datum_depth = missing;
    woc = NaN,
    goc = NaN,
    wgc = NaN,
    pc_woc = 0.0,
    pc_goc = 0.0,
    pc_wgc = 0.0,
    temperature = missing,
    temperature_vs_depth = missing,
    composition = missing,
    composition_vs_depth = missing,
    liquid_composition = missing,
    liquid_composition_vs_depth = missing,
    vapor_composition = missing,
    vapor_composition_vs_depth = missing,
    density_function = missing,
    rs = 0.0,
    rs_vs_depth = missing,
    rv = 0.0,
    rv_vs_depth = missing,
    pvtnum = 1,
    satnum = 1,
    cells = missing,
    kwarg...
)
```


Set uip equilibriation region for a reservoir model. The region is defined by the datum pressure and depth, and the water-oil, gas-oil, and water-gas contacts. The contacts will be used to determine the phase distribution and initial pressure in the region. The region can be further specified by temperature, composition, density, and rs/rv functions. Most entries can either be specified as a function of depth or as a constant value. Additional keyword arguments are passed onto the `equilibriate_state` function.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/state0.jl#L21-L56" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.equilibriate_state' href='#JutulDarcy.equilibriate_state'><span class="jlbinding">JutulDarcy.equilibriate_state</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
equilibriate_state(model, contacts)
```


Equilibrates the state of the given model based on the provided contacts.

**Arguments**
- `model`: The model whose state needs to be equilibrated.
  
- `contacts`: The nph contact depths.
  

**Keyword Arguments**
- `datum_depth`: The reference depth for the datum.
  
- `datum_pressure`: The pressure at the datum depth.
  
- `cells`: The cells to be equilibrated.
  
- `rs`: Solution gas-oil ratio (blackoil).
  
- `rv`: Vapor-oil ratio (blackoil).
  
- `composition`: The composition vs depth (compositional).
  
- `kwarg`: Additional keyword arguments.
  

**Returns**
- The equilibrated state of the model.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/init/init.jl#L2-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Simulation {#Simulation}

Simulating is done by either setting up a reservoir simulator and then simulating, or by using the convenience function that automatically sets up a simulator for you.

There are a number of different options available to tweak the tolerances, timestepping and behavior of the simulation. It is advised to read through the documentation in detail before running very large simulations.
<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.simulate_reservoir' href='#JutulDarcy.simulate_reservoir'><span class="jlbinding">JutulDarcy.simulate_reservoir</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
simulate_reservoir(state0, model, dt;
    parameters = setup_parameters(model),
    restart = false,
    forces = setup_forces(model),
    kwarg...
)
simulate_reservoir(case;
    kwarg...
)
```


Convenience function for simulating a reservoir model. This function internally calls [`setup_reservoir_simulator`](/man/highlevel#JutulDarcy.setup_reservoir_simulator), simulates the problem and returns a [`ReservoirSimResult`](/man/highlevel#JutulDarcy.ReservoirSimResult). Keyword arguments are passed onto [`setup_reservoir_simulator`](/man/highlevel#JutulDarcy.setup_reservoir_simulator) and are documented in that function.

You can optionally unpack this result into the most typical desired outputs:

`wellsols, states = simulate_reservoir(...)`

where `wellsols` contains the well results and `states` the reservoir results (pressure, saturations and so on, in each cell of the reservoir domain).

**Examples**

You can restart/resume simulations by both providing the `output_path` argument and the `restart` argument:

```
# Automatically restart from last solved step and returning the outputs if the simulation was already solved.
result = simulate_reservoir(state0, model, dt, output_path = "/some/path", restart = true)

# Restart from step 5
result = simulate_reservoir(state0, model, dt, output_path = "/some/path", restart = 5)

# Start from the beginning (default)
result = simulate_reservoir(state0, model, dt, output_path = "/some/path", restart = false)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L1058-L1095" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.setup_reservoir_simulator' href='#JutulDarcy.setup_reservoir_simulator'><span class="jlbinding">JutulDarcy.setup_reservoir_simulator</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
setup_reservoir_simulator(models, initializer, parameters = nothing; <keyword arguments>)
```


**Arguments**
- `models`: either a single model or a Dict with the key :Reservoir for multimodels
  
- `initializer`: used to setup state0, must be compatible with `model`
  
- `parameters`: initialized parameters, must be compatible with `model` if provided
  

**Keyword arguments**
- `split_wells`: Add facility model to each well (needed for domain decomposition and MPI solves)
  
- `assemble_wells_together`: Option to split wells into multiple sparse matrices (false argument experimental)
  
- `specialize=false`: use deep specialization of storage for faster execution, but significantly more compile time
  

Additional keyword arguments are documented in the version of this function that uses `JutulCase` as the input.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L756-L770" target="_blank" rel="noreferrer">source</a></Badge>



```julia
setup_reservoir_simulator(case::JutulCase; <keyword arguments>)
```


**Keyword arguments**
- `mode=:default`: Mode used for solving. Can be set to `:mpi` if running in MPI mode together with HYPRE, PartitionedArrays and MPI in your environment.
  
- `method=:newton`: Can be `:newton`, `:nldd` or `:aspen`. Newton is the most tested approach and `:nldd` can speed up difficult models. The `:nldd` option enables a host of additional options (look at the simulator config for more details).
  
- `presolve_wells=false`: Solve wells before each linearization. Can improve convergence for models with multisegment wells.
  

**Linear solver options**
- `linear_solver=:bicgstab`: iterative solver to use (provided model supports it). Typical options are `:bicgstab` or `:gmres` Can alternatively pass a linear solver instance.
  
- `precond=:cpr`: preconditioner for iterative solver. For larger problems, CPR variants are recommended. In order of strength and cost:
  - `:cpr` standard Constrained-Pressure-Residual with ILU(0) second stage (strong preconditioner)
    
  - `:cprw` CPRW with ILU(0) second stage. Faster for problems with wells (strong preconditioner)
    
  - `:ilu0` block-incomplete-LU (intermediate strength preconditioner)
    
  - `:spai0`: Sparse Approximate Inverse of lowest order (weak preconditioner)
    
  - `jacobi`: Jacobi preconditioner (weak preconditioner)
    
  
- `rtol=nothing`: relative tolerance for linear solver. If set to `nothing`, the default tolerance for the preconditioner is used, which is 5e-3 for CPR variants and 1e-2 for smoothers.
  
- `linear_solver_arg`: `Dict` containing additional linear solver arguments.
  

**Timestepping options**
- `initial_dt=si_unit(:day)`: initial timestep in seconds (one day by default)
  
- `target_ds=Inf`: target saturation change over a timestep used by timestepper.
  
- `target_dz=Inf`: target mole fraction change over a timestep used by timestepper (compositional only).
  
- `target_its=8`: target number of nonlinear iterations per time step
  
- `offset_its=1`: dampening parameter for time step selector where larger values lead to more pessimistic estimates.
  
- `timesteps=:auto`: Set to `:auto` to use automatic timestepping, `:none` for no automatic timestepping (i.e. try to solve exact report steps)
  
- `max_timestep=si_unit(:year)`: Maximum internal timestep used in solver.
  
- `min_timestep=0.0`: Minimum internal timestep used in solver.
  

**Convergence criterions**
- `tol_cnv=1e-3`: maximum allowable point-wise error (volume-balance)
  
- `tol_mb=1e-7`: maximum alllowable integrated error (mass-balance)
  
- `tol_cnv_well=10*tol_cnv`: maximum allowable point-wise error for well node (volume-balance)
  
- `tol_mb_well=1e4*tol_mb`: maximum alllowable integrated error for well node (mass-balance)
  
- `inc_tol_dp_abs=Inf`: Maximum allowable pressure change (absolute)
  
- `inc_tol_dp_rel=Inf`: Maximum allowable pressure change (absolute)
  
- `inc_tol_dz=Inf`: Maximum allowable composition change (compositional only).
  

**Inherited keyword arguments**

Additional keyword arguments come from the base Jutul simulation framework. We list a few of the most relevant entries here for convenience:
- `info_level = 0`: Output level. Set to 0 for minimal output, -1 for no output and 1 or more for increasing verbosity.
  
- `output_path`: Path to write output to.
  
- `max_nonlinear_iterations=15`: Maximum Newton iterations before a time-step is cut.
  
- `min_nonlinear_iterations=1`: Minimum number of Newtons to perform before checking convergence.
  
- `relaxation=Jutul.NoRelaxation()`: Dampening used for solves. Can be set to `Jutul.SimpleRelaxation()` for difficult models. Equivialent option is to set `true` for relaxation and `false` for no relaxation.
  
- `failure_cuts_timestep=true`: Cut timestep instead of throwing an error when numerical issues are encountered (e.g. linear solver divergence).
  
- `max_timestep_cuts=25`: Maximum number of timestep cuts before a solver gives up. Note that when using dynamic timestepping, this in practice defines a minimal timestep, with more than the prescribed number of cuts being allowed if the timestep is dynamically increased after cutting.
  
- `timestep_max_increase=10.0`: Max allowable factor to increase time-step by. Overrides any choices made in dynamic step selection.
  
- `timestep_max_decrease=0.1`: Max allowable factor to decrease time-step by. Overrides any choices made in dynamic step selection.
  
- `tol_factor_final_iteration=1.0`: If set to a value larger than 1.0, the final convergence check before a time-step is cut is relaxed by multiplying all tolerances with this value. Warning: Setting it to a large value can have severe impact on numerical accuracy. A value of 1 to 10 is typically safe if your default tolerances are strict.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/utils.jl#L812-L900" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='JutulDarcy.ReservoirSimResult' href='#JutulDarcy.ReservoirSimResult'><span class="jlbinding">JutulDarcy.ReservoirSimResult</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ReservoirSimResult(model, result::Jutul.SimResult, forces, extra = Dict(); kwarg...)
```


Create a specific reservoir simulation results that contains well curves, reservoir states, and so on. This is the return type from `simulate_reservoir`.

A `ReservoirSimResult` can be unpacked into well solutions, reservoir states and reporting times:

```julia
res_result::ReservoirSimResult
ws, states, t = res_result
```


**Fields**
- `wells`
  
- `states`
  
- `time`
  
- `result`
  
- `extra`
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/sintefmath/JutulDarcy.jl/blob/d6c48cb7b990b0808b12f3f1a008f1429a0152ed/src/types.jl#L639-L656" target="_blank" rel="noreferrer">source</a></Badge>

</details>

