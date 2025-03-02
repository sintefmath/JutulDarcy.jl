# # Gravity circulation with CPR preconditioner
# This example demonstrates a more complex gravity driven instability. The
# problem is a bit larger than the [Gravity segregation example](@ref), and is
# therefore set up using the high level API that automatically sets up an
# iterative linear solver with a constrained pressure residual (CPR)
# preconditioner and automatic timestepping.
#
# The high level API uses the more low level `Jutul` API seen in the other
# examples under the hood and makes more complex problems easy to set up. The
# same data structures and functions are used, allowing for deep customization
# if the defaults are not appropriate.
using JutulDarcy
using Jutul
using GLMakie
cmap = :seismic
nx = nz = 100;
# ## Define the domain
D = 10.0
g = CartesianMesh((nx, 1, nz), (D, 1.0, D))
domain = reservoir_domain(g);
# ## Set up model and properties
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)
p0 = 100*bar
rhoLS = 1000.0*kg/meter^3 # Definition of fluid phases
rhoVS = 500.0*kg/meter^3
cl, cv = 1e-5/bar, 1e-4/bar
L, V = LiquidPhase(), VaporPhase()
sys = ImmiscibleSystem([L, V])
model, parameters = setup_reservoir_model(domain, sys)
density = ConstantCompressibilityDensities(sys, p0, [rhoLS, rhoVS], [cl, cv]) # Replace density with a lighter pair
replace_variables!(model, PhaseMassDensities = density);
kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 3.0])
replace_variables!(model, RelativePermeabilities = kr)

# ### Define initial saturation
# Set the left part of the domain to be filled by the vapor phase and the heavy
# liquid phase in the remainder. To do this, we grab the cell centroids in the x
# direction from the domain, reshape them to the structured mesh we are working
# on and define the liquid saturation from there.
c = domain[:cell_centroids]
x = reshape(c[1, :], nx, nz)

sL = zeros(nx, nz)
plane = D/2.0
for i in 1:nx
    for j = 1:nz
        X = x[i, j]
        sL[i, j] = clamp(Float64(X > plane), 0, 1)
    end
end
heatmap(sL, colormap = cmap, axis = (title = "Initial saturation",))
# ### Set up initial state
sL = vec(sL)'
sV = 1 .- sL
s0 = vcat(sV, sL)
state0 = setup_reservoir_state(model, Pressure = p0, Saturations = s0)

# ### Set the viscosity of the phases
# By default, viscosity is a parameter and can be set per-phase and per cell.
μ = parameters[:Reservoir][:PhaseViscosities]
@. μ[1, :] = 1e-3
@. μ[2, :] = 5e-3
# Convert time-steps from days to seconds
timesteps = repeat([10.0*3600*24], 20)
_, states, = simulate_reservoir(state0, model, timesteps, parameters = parameters);
# ## Plot results
# Plot initial saturation
tmp = reshape(state0[:Reservoir][:Saturations][1, :], nx, nz)
f = Figure()
ax = Axis(f[1, 1], title = "Before")
heatmap!(ax, tmp, colormap = cmap)
# Plot intermediate saturation
tmp = reshape(states[length(states) ÷ 2][:Saturations][1, :], nx, nz)
ax = Axis(f[1, 2], title = "Half way")
hm = heatmap!(ax, tmp, colormap = cmap)
# Plot final saturation
tmp = reshape(states[end][:Saturations][1, :], nx, nz)
ax = Axis(f[1, 3], title = "After")
hm = heatmap!(ax, tmp, colormap = cmap)
Colorbar(f[1, 4], hm)
f
