# # Gravity segregation example
# The simplest type of porous media simulation problem to set up that is not
# trivial is the transition to equilibrium from an unstable initial condition.
# Placing a heavy fluid on top of a lighter fluid will lead to the heavy fluid
# moving down while the lighter fluid moves up. 
#
# ## Problem set up
# We define a simple 1D gravity column with an approximate 10-1 ratio in density
# between the two compressible phases and let it simulate until equilibrium is
# reached. We begin by defining the reservoir domain itself.
using JutulDarcy, Jutul
nc = 100
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)

g = CartesianMesh((1, 1, nc), (1.0, 1.0, 10.0))
domain = reservoir_domain(g, permeability = 1.0*Darcy);
# ## Fluid properties
# Define two phases liquid and vapor with a 10-1 ratio reference densities and
# set up the simulation model.
p0 = 100*bar

rhoLS = 1000.0*kg/meter^3
rhoVS = 100.0*kg/meter^3
cl, cv = 1e-5/bar, 1e-4/bar
L, V = LiquidPhase(), VaporPhase()
sys = ImmiscibleSystem([L, V])
model = SimulationModel(domain, sys);
# ### Definition for phase mass densities
# Replace default density with a constant compressibility function that uses the
# reference values at the initial pressure.
density = ConstantCompressibilityDensities(sys, p0, [rhoLS, rhoVS], [cl, cv])
set_secondary_variables!(model, PhaseMassDensities = density);
# ### Set up initial state
# Put heavy phase on top and light phase on bottom. Saturations have one value
# per phase, per cell and consequently a per-cell instantiation will require a
# two by number of cells matrix as input. We also set up time-steps for the
# simulation, using the provided conversion factor to convert days into seconds.
nl = nc ÷ 2
sL = vcat(ones(nl), zeros(nc - nl))'
s0 = vcat(sL, 1 .- sL)
state0 = setup_state(model, Pressure = p0, Saturations = s0)
timesteps = repeat([0.02]*day, 150);
# ## Perform simulation
# We simulate the system using the default linear solver and otherwise default
# options. Using `simulate` with the default options means that no dynamic
# timestepping will be used, and the simulation will report on the exact 150
# steps defined above.
states, report = simulate(state0, model, timesteps);
# ## Plot results
# Plot the saturations of the liquid phase at three different timesteps: The
# initial, unstable state, an intermediate state where fluid exchange between
# the top and bottom is initiated, and the final equilibrium state where the
# phases have swapped places.
using GLMakie
fig = Figure()
function plot_sat!(ax, state)
    plot_cell_data!(ax, g, state[:Saturations][1, :],
        colorrange = (0.0, 1.0),
        colormap = :seismic
    )
end
ax1 = Axis3(fig[1, 1], title = "Initial state", aspect = (1, 1, 4.0))
plot_sat!(ax1, state0)

ax2 = Axis3(fig[1, 2], title = "Intermediate state", aspect = (1, 1, 4.0))
plot_sat!(ax2, states[25])

ax3 = Axis3(fig[1, 3], title = "Final state", aspect = (1, 1, 4.0))
plt = plot_sat!(ax3, states[end])
Colorbar(fig[1, 4], plt)
fig
# ### Plot time series
# The 1D nature of the problem allows us to plot all timesteps simultaneously in
# 2D. We see that the heavy fluid, colored blue, is initially at the top of the
# domain and the lighter fluid is at the bottom. These gradually switch places
# until all the heavy fluid is at the lower part of the column.
tmp = vcat(map((x) -> x[:Saturations][1, :]', states)...)
f = Figure()
ax = Axis(f[1, 1], xlabel = "Time", ylabel = "Depth", title = "Gravity segregation")
hm = heatmap!(ax, tmp, colormap = :seismic)
Colorbar(f[1, 2], hm)
f
