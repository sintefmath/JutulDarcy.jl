# # Introduction to wells

# This example demonstrates how to set up a 3D domain with a layered
# permeability field, define wells and solve a simple production-injection
# schedule. We begin by loading the `Jutul` package that contains generic
# features like grids and linear solvers and the `JutulDarcy` package itself.

# ## Preliminaries
using JutulDarcy, Jutul

# `JutulDarcy` uses SI units internally. It is therefore convenient to define a
# few constants at the start of the script to have more managable numbers later
# on.
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day);

# ## Defining a porous medium
# We start by defining the static part of our simulation problem -- the porous medium itself.
# ### Defining the grid

# The first step is to create a grid for our simulation domain. We make a tiny 5
# by 5 grid with 4 layers that discretizes a physical domain of 2000 by 1500 by
# 50 meters.
nx = ny = 5
nz = 4
dims = (nx, ny, nz)
g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
#-
# ### Adding properties and making a domain
# The grid by itself does not fully specify a porous medium. For that we need to
# specify the permeability in each cell and the porosity. Permeability, often
# denoted by a positive-definite tensor K, describes the relationship between a
# pressure gradient and the flow through the medium. Porosity is a dimensionless
# number between 0 and 1 that describes how much of the porous medium is void
# space where fluids can be present. The assumption of Darcy flow becomes less
# reasonable for high porosity values and the flow equations break down at zero
# porosity. A porosity of 0.2 is then a safe choice.

# Jutul uses the `DataDomain` type to store a domain/grid together with data.
# For porous media simulation, `JutulDarcy` includes a convenience function
# `reservoir_domain` that contains defaults for permeability and porosity. We
# specify the permeability per-cell with varying values per layer in the
# vertical direction and a single porosity value for all cells that the function
# will expand for us. From the output, we can see that basic geometry primitives
# are also automatically added:
nlayer = nx*ny # Cells in each layer
K = vcat(
    fill(0.65, nlayer),
    fill(0.3, nlayer),
    fill(0.5, nlayer),
    fill(0.2, nlayer)
    )*Darcy

domain = reservoir_domain(g, permeability = K, porosity = 0.2)

# ## Defining wells
# Now that we have a porous medium with all static properties set up, it is time
# to introduce some driving forces. Jutul assumes no-flow boundary conditions on
# all boundary faces unless otherwise specified so we can go ahead and add wells
# to the model. 

# ### A vertical producer well
# We will define two wells: A first well is named "Producer" and is a vertical
# well positioned at `(1, 1)`. By default, the [`setup_vertical_well`](@ref)
# function perforates all layers in the model.
Prod = setup_vertical_well(domain, 1, 1, name = :Producer);

# ### A single-perforation injector
# We also define an injector by [`setup_well`](@ref). This function allows us to
# pass a vector of either cell indices or tuples of logical indices that the
# well trajectory will follow. We setup the injector in the upper left corner.

Inj = setup_well(domain, [(nx, ny, 1)], name = :Injector);
# ## Choosing a fluid system
# To solve multiphase flow with our little toy reservoir we need to pick a fluid
# system. The type of system determines what physical effects are modelled, what
# parameters are required and the runtime and accuracy of the resulting
# simulation. The choice is in practice a trade-off between accuracy, runtime
# and available data that should be informed by modelling objectives. In this
# example our goal is to understand how to set up a simple well problem and the
# [`ImmiscibleSystem`](@ref) requires a minimal amount of input. We define
# liquid and gas phases and their densities at some reference conditions and
# instantiate the system.
## Set up a two-phase immiscible system and define a density secondary variable
phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0
rhoGS = 100.0
rhoS = [rhoLS, rhoGS] .* kg/meter^3
sys = ImmiscibleSystem(phases, reference_densities = rhoS)
# ### Creating the model
# The same fluid system can be used for both flow inside the wells and the
# reservoir. JutulDarcy treats wells as first-class citizens and models flow
# inside the well bore using the same fluid description as the reservoir, with
# modified equations to account for the non-Darcy velocities. We call the
# utility function that sets up all of this for us:
model, parameters = setup_reservoir_model(domain, sys, wells = [Inj, Prod])
model
# The model is an instance of the [`MultiModel`](@ref) from `Jutul` where a
# submodel is defined for the reservoir, each of the wells and the facility that
# controls both wells. In addition we can see the cross-terms that couple these
# wells together. If we want to see more details on how either of these are set
# up, we can display for example the reservoir model.
reservoir = model[:Reservoir]

# We can see that the model contains primary variables, secondary variables
# (sometimes referred to as properties) and static parameters in addition to the
# system we already set up. These can be replaced or modified to alter the
# behavior of the system.

# ### Replace the density function with our custom version
# Let us change the definition of phase mass densities for our model. We'd like
# to model our liquid phase as weakly compressible and the vapor phase with more
# significant compressibility. A common approach is to define densities
# ``\rho_\alpha^s`` at some reference pressure ``p_r`` and use a phase
# compressibility ``c_\alpha`` to extrapolate around that known value.
#
# ``\rho_\alpha (p) = \rho_\alpha^s \exp((p - p_r)c_\alpha)``
#
# This is already implement in Jutul and we simply need to instantiate the variable definition:
c = [1e-6/bar, 1e-4/bar]
ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
# Before replacing it in the model. This change will propagate to all submodels
# that have a definition given for PhaseMassDensities, including the wells. The
# facility, which does not know about densities, will ignore it.
replace_variables!(model, PhaseMassDensities = ρ);
# This concludes the setup of the model.
# ## Set up initial state
# The model is time-dependent and requires initial conditions. For the
# immiscible model it is sufficient to specify the reference phase pressure and
# the saturations for both phases, summed up to one. These can be specified per
# cell or one for the entire grid. Specifying a single pressure for the entire
# model is not very realistic, but should be fine for our simple example. The
# initial conditions will equilibrate themselves from gravity fairly quickly.
state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])

# ## Set up report time steps and injection rate
# We create a set of time-steps. These are report steps where the solution will
# be reported, but the simulator itself will do internal subdivision of time
# steps if these values are too coarse for the solvers. We also define an
# injection rate of a full pore-volume (at reference conditions) of gas.
dt = repeat([30.0]*day, 12*5)
pv = pore_volume(model, parameters)
inj_rate = sum(pv)/sum(dt)
# ## Set up well controls
# We then set up a total rate target (positive value for injection) together
# with a corresponding injection control that specifies the mass fractions of
# the two components/phases for pure gas injection, with surface density given
# by the known gas density. The producer operates at a fixed bottom hole
# pressure. These are given as a `Dict` with keys that correspond to the well
# names.
rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
bhp_target = BottomHolePressureTarget(50*bar)
P_ctrl = ProducerControl(bhp_target)
controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl
## Set up the forces
# Set up forces for the whole model. For this example, all other forces than the
# well controls are defaulted (amounting to no-flow for the reservoir). Jutul
# supports either a single set of forces for the entire simulation, or a vector
# of equal length to `dt` with varying forces. Reasonable operational limits for
# wells are also set up by default.
forces = setup_reservoir_forces(model, control = controls)

# ## Simulate the model
# We are finally ready to simulate the model for the given initial state
# `state0`, report steps `dt`, `parameters` and forces. As the model is small,
# barring any compilation time, this should run in about 300 ms.
result = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces);
# ### Unpacking the result
# The result contains a lot of data. This can be unpacked to get the most
# typical desired outputs: Well responses, reservoir states and the time they
# correspond to. 
wd, states, t = result;
# We could in fact equally well have written
# `wd, states, t = simulate_reservoir(...)`
# to arrive at the same result.

# ## Plot the producer responses
# We load a plotting package to plot the wells.
using GLMakie
# ## Plot the surface rates at the producer
# We observe that the total rate does not vary much, but the composition changes
# from liquid to gas as the front propagate through the domain and hits the
# producer well.
# Gas rates:
qg = wd[:Producer][:grat];
# Total rate:
qt = wd[:Producer][:rate];
# Compute liquid rate and plot:
ql = qt - qg
x = t/day
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time (days)",
    ylabel = "Rate (m³/day)",
    title = "Well production rates"
)
lines!(ax, x, abs.(qg).*day, label = "Gas")
lines!(ax, x, abs.(ql).*day, label = "Liquid")
lines!(ax, x, abs.(qt).*day, label = "Total")
axislegend(position = :rb)
fig
# ## Plot bottom hole pressure of the injector
# The pressure builds during injection, until the gas breaks through to the
# other well.
bh = wd[:Injector][:bhp]
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time (days)",
    ylabel = "Bottom hole pressure (bar)",
    title = "Injector bottom hole pressure"
)
lines!(ax, x, bh./bar)
fig
# ## Plot the well results in the interactive viewer
# Note that this will open a new window with the plot.
plot_well_results(wd, resolution = (800, 500))
# ## Plot the reservoir and final gas saturation field
plot_cell_data(g, states[end][:Saturations][1, :], colormap = :seismic)
