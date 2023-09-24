# # Gravity segregation example
# The simplest type of porous media simulation problem to set up that is not
# trivial is the transition to equilibrium from an unstable initial condition.
# Placing a heavy fluid on top of a lighter fluid will lead to the heavy fluid
# moving down while the lighter fluid moves up. 
#
# ## Problem set up
# We define a simple 1D gravity column with an approximate 10-1 ratio in density
# between the two compressible phases and let it simulate until equilibrium is
# reached.
using JutulDarcy, Jutul
nc = 100
domain = get_1d_reservoir(nc, z_max = 1)
#-
# ## Fluid properties
# Define two phases liquid and vapor with a 10-1 ratio reference densities and
# set up the simulation model.
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)
p0 = 100*bar

rhoLS = 1000.0*kg/meter^3
rhoVS = 100.0*kg/meter^3
cl, cv = 1e-5/bar, 1e-4/bar
L, V = LiquidPhase(), VaporPhase()
sys = ImmiscibleSystem([L, V])
model = SimulationModel(domain, sys)
#-
# ### Definition for phase mass densities
# Replace default density with a constant compressibility function that uses the
# reference values at the initial pressure.
density = ConstantCompressibilityDensities(sys, p0, [rhoLS, rhoVS], [cl, cv])
set_secondary_variables!(model, PhaseMassDensities = density)
# ### Set up initial state
# Put heavy phase on top and light phase on bottom. Saturations have one value
# per phase, per cell and consequently a per-cell instantiation will require a
# two by number of cells matrix as input.
nl = nc ÷ 2
sL = vcat(ones(nl), zeros(nc - nl))'
s0 = vcat(sL, 1 .- sL)
state0 = setup_state(model, Pressure = p0, Saturations = s0)
# Convert time-steps from days to seconds
timesteps = repeat([0.02]*day, 150)
#-
## Perform simulation
states, report = simulate(state0, model, timesteps, info_level = -1)
# ## Plot results
# The 1D nature of the problem allows us to plot all timesteps simultaneously in
# 2D. We see that the heavy fluid, colored blue, is initially at the top of the
# domain and the lighter fluid is at the bottom. These gradually switch places
# until all the heavy fluid is at the lower part of the column.
using CairoMakie
tmp = vcat(map((x) -> x[:Saturations][1, :]', states)...)
f = Figure()
ax = Axis(f[1, 1], xlabel = "Time", ylabel = "Depth", title = "Gravity segregation")
hm = heatmap!(ax, tmp, colormap = :seismic)
Colorbar(f[1, 2], hm)
f
