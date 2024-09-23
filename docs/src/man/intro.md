# Getting started

## Installing JutulDarcy

To get started with JutulDarcy, you must install [the Julia programming language](https://julialang.org/). We recommend the latest stable release, but at least version 1.9 is required. Julia uses environments to manage packages. If you are not familiar with this concept, we recommend the [Pkg documentation](https://pkgdocs.julialang.org/v1/environments/).

### Setting up an environment

To set up an environment you can create a folder and open Julia with that project as a runtime argument. The environment is persistent: If you start Julia in the same folder with the same project argument, the same packages will be installed. For this reason, it is sufficient to do this process once.

```bash
mkdir jutuldarcy_env
cd jutuldarcy_env
julia --project=.
```

```julia
using Pkg
Pkg.add("JutulDarcy")
```

You can then run any of the examples in the [`examples`](https://github.com/sintefmath/JutulDarcy.jl/tree/main/examples) directory by including them. The examples are sorted by complexity. We suggest you start with [Gravity segregation example](Intro: Gravity segregation example).

To generate a folder that contains the examples locally, you can run the following code to create a folder `jutuldarcy_examples` in your current working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples()
```

Alternatively, a folder can be specified if you want the examples to be placed outside your present working directory:

```julia
using JutulDarcy
generate_jutuldarcy_examples("/home/username/")
```

### Adding additional packages

We also rely heavily on the Jutul base package for some functionality, so we recommend that you also install it together with JutulDarcy:

```julia
Pkg.add("Jutul")
```

If you want the plotting used in the examples, you need this:

```julia
Pkg.add("GLMakie") # 3D and interactive visualization
```

Some examples and functionaltiy also make use of additional packages:

```julia
Pkg.add("Optim") # Optimization library
Pkg.add("HYPRE") # Better linear solver
Pkg.add("GeoEnergyIO") # Parsing input files
```

## A short motivational example with two wells

After a bit of time has passed compiling the packages, you are now ready to use JutulDarcy. There are a number of examples included in this manual, but we include a brief example here that briefly demonstrates the key concepts.

### Setting up the domain

We set up a simple Cartesian Mesh that is converted into a reservoir domain with permeability and porosity values. We then use this domain to set up two wells: One vertical well for injection and a single perforation producer:

````@example intro_ex
using JutulDarcy, Jutul
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)
nx = ny = 25
nz = 10
cart_dims = (nx, ny, nz)
physical_dims = (1000.0, 1000.0, 100.0).*meter
g = CartesianMesh(cart_dims, physical_dims)
domain = reservoir_domain(g, permeability = 0.3Darcy, porosity = 0.2)
Injector = setup_vertical_well(domain, 1, 1, name = :Injector)
Producer = setup_well(domain, (nx, ny, 1), name = :Producer)
# Show the properties in the domain
domain
````

### Setting up a fluid system

We select a two-phase immicible system by declaring that the liquid and vapor phases are present in the model. These are assumed to have fixed densitys of 1000 and 100 kilograms per meters cubed at some reference pressure and temperature conditions.

````@example intro_ex
phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0kg/meter^3
rhoGS = 100.0kg/meter^3
reference_densities = [rhoLS, rhoGS]
sys = ImmiscibleSystem(phases, reference_densities = reference_densities)
````

### Setting up the model

We now have everything we need to set up a model. We call the setup function and get out the model together with the parameters - numerical input values that are static throughout the simulation. These are automatically computed from the domain's geometry, permeability and porosity.

````@example intro_ex
model, parameters = setup_reservoir_model(domain, sys, wells = [Injector, Producer])
model
````

The model has a set of default secondary variables (properties) that are used to compute the flow equations. We can have a look at the reservoir model to see what the defaults are for the Darcy flow part of the domain:

````@example intro_ex
reservoir_model(model)
````

The secondary variables can be swapped out, replaced and new variables can be added with arbitrary functional dependencies thanks to Jutul's flexible setup for automatic differentiation. Let us adjust the defaults by replacing the relative permeabilities with Brooks-Corey functions and the phase density functions by constant compressibilities:

````@example intro_ex
c = [1e-6, 1e-4]/bar
density = ConstantCompressibilityDensities(
    p_ref = 100*bar,
    density_ref = reference_densities,
    compressibility = c
)
kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 3.0])
replace_variables!(model, PhaseMassDensities = density, PhaseRelativePermeability = kr)
nothing #hide
````

### Initial state

Now that we are happy with our model setup, we can designate an initial state. For simplicity, we set up a constant pressure reservoir filled with the liquid phase.

````@example intro_ex
state0 = setup_reservoir_state(model,
    Pressure = 120bar,
    Saturations = [1.0, 0.0]
)
````

### Setting up timesteps and well controls

We set up reporting timesteps. These are the intervals that the simulator gives out outputs. The simulator may use shorter steps internally, but will always hit these points in the output.

````@example intro_ex
nstep = 25
dt = fill(365.0day, nstep)
nothing #hide
````

We next set up a rate target with a high amount of gas injected into the model. This is not fully realistic, but will give some nice and dramatic plots for our example later on.

````@example intro_ex
pv = pore_volume(model, parameters)
inj_rate = 1.5*sum(pv)/sum(dt)
rate_target = TotalRateTarget(inj_rate)
````

The producer is set to operate at a fixed pressure:

````@example intro_ex
bhp_target = BottomHolePressureTarget(100bar)
````

We can finally set up forces for the model. Note that while JutulDarcy supports time-dependent forces and limits for the wells, we keep this example as simple as possible.

````@example intro_ex
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
P_ctrl = ProducerControl(bhp_target)
controls = Dict(:Injector => I_ctrl, :Producer => P_ctrl)
forces = setup_reservoir_forces(model, control = controls)
nothing #hide
````

### Simulate and analyze results

We call the simulation with our initial state, our model, the timesteps, the forces and the parameters:

````@example intro_ex
wd, states, t = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces)
````

We can interactively look at the well results in the command line:

````@example intro_ex
wd(:Producer)
````

Let us look at the pressure evolution in the injector:

````@example intro_ex
wd(:Injector, :bhp)
````

If we have a plotting package available, we can visualize the results too:

````@example intro_ex
using GLMakie
grat = wd[:Producer, :grat]
lrat = wd[:Producer, :lrat]
bhp = wd[:Injector, :bhp]
fig = Figure(size = (1200, 400))
ax = Axis(fig[1, 1],
    title = "Injector",
    xlabel = "Time / days",
    ylabel = "Bottom hole pressure / bar")
lines!(ax, t/day, bhp./bar)
ax = Axis(fig[1, 2],
    title = "Producer",
    xlabel = "Time / days",
    ylabel = "Production rate / m³/day")
lines!(ax, t/day, abs.(grat).*day)
lines!(ax, t/day, abs.(lrat).*day)
fig
````

Interactive visualization of the 3D results is also possible if GLMakie is loaded:

````@example intro_ex
plot_reservoir(model, states, key = :Saturations, step = 3)
````
