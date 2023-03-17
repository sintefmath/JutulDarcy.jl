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
function solve_gravity_column(nc = 100, tstep = repeat([0.02], 150))
    domain = get_1d_reservoir(nc, z_max = 1)
    nc = number_of_cells(domain)
    bar = 1e5
    p0 = 100*bar
    rhoLS, rhoVS = 1000.0, 100.0 # Definition of fluid phases
    cl, cv = 1e-5/bar, 1e-4/bar
    L, V = LiquidPhase(), VaporPhase()
    sys = ImmiscibleSystem([L, V])
    model = SimulationModel(domain, sys)
    density = ConstantCompressibilityDensities(sys, p0, [rhoLS, rhoVS], [cl, cv]) # Replace density with a lighter pair
    set_secondary_variables!(model, PhaseMassDensities = density)
    nl = nc ÷ 2
    sL = vcat(ones(nl), zeros(nc - nl))'
    s0 = vcat(sL, 1 .- sL) # Put heavy phase on top and light phase on bottom
    state0 = setup_state(model, Pressure = p0, Saturations = s0)
    timesteps = tstep*3600*24 # Convert time-steps from days to seconds
    states, report = simulate(state0, model, timesteps, info_level = -1)
    return states, model, report
end

## Perform simulation
states, model, report = solve_gravity_column();
nothing

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
