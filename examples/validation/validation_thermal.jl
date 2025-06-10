# # Aquifer thermal energy storage (ATES) validation
# This example validates JutulDarcy's thermal solver against results from a
# commercial simulator. The test case is a simple ATES model with a single pair
# of wells (hot / cold) where the cold well is used for pressure support with a
# mirrored injection rate. Towards the later part of the schedule, the cold well
# reinjects water that is a higher temperature than the background.
#
# This case is completely specified in the `ATES_TEST.DATA` file which was
# provided by TNO. The model is a structured mesh with 472 500 active cells.
using Jutul, JutulDarcy, GeoEnergyIO, DelimitedFiles, HYPRE, GLMakie
basepth = GeoEnergyIO.test_input_file_path("ATES_TEST")
data = parse_data_file(joinpath(basepth, "ATES_TEST.DATA"))

wdata, wheader = readdlm(joinpath(basepth, "wells.txt"), ',', header = true)
cdata, cheader = readdlm(joinpath(basepth, "cells.txt"), ',', header = true)
case = setup_case_from_data_file(data)
ws, states, t_seconds = simulate_reservoir(case, info_level = 1);
# ## Plot the reservoir and monitor points
# We will monitor points close to the warm and cold wells for comparison with a
# commercial simulator. These can be identified by their IJK triplets, and we
# plot these in orange and blue.
reservoir = reservoir_domain(case)
G = physical_representation(reservoir)

cell_warm1 = cell_index(G, (105,75,6))
cell_warm2 = cell_index(G, (106,75,6))

cell_cold1 = cell_index(G, (50,75,6))
cell_cold2 = cell_index(G, (51,75,6))

fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1], zreversed = true, azimuth = 4.55, elevation = 0.2)
for (wname, w) in get_model_wells(case)
    plot_well!(ax, G, w)
end
Jutul.plot_mesh_edges!(ax, G, alpha = 0.1)
plot_mesh!(ax, G, color = :orange, cells = [cell_warm1, cell_warm2])
plot_mesh!(ax, G, color = :blue, cells = [cell_cold1, cell_cold2])
fig
# ## Plot the permeability
# The model contains a high permeable aquifer layer in the middle of the model.
# The high permeable layer conducts the flow between the wells, and the low
# permeable layers above and below the aquifer layer conduct heat.
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1], zreversed = true, azimuth = 4.55, elevation = 0.2)
for (wname, w) in get_model_wells(case)
    plot_well!(ax, G, w)
end
plt = plot_cell_data!(ax, G, reservoir[:permeability][1, :]./si_unit(:darcy),
    shading = NoShading,
    colormap = :thermal
)
Colorbar(fig[2, 1], plt, label = "Horizontal permeability (darcy)", vertical = false)
fig
# ## Plot the water rate in the wells
# The water rate in the wells is shown below. The well rates are mirrored in
# that the cold well injects the same amount of water as the warm well produces
# and vice versa. Initially, the warm well injects water and the cold well
# produces to maintain aquifer pressure. Later on, the warm well produces warm
# water and the cold well injects utilized cold water at a slightly higher
# temperature than that of the reservoir to balance the pressure.
day = si_unit(:day)
t_jutul = t_seconds./day
wrat_cold = ws[:COLD][:wrat]*day
wrat_warm = ws[:WARM][:wrat]*day
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time elapsed (days)", ylabel = "Water rate (m³/day)", title = "Water rate in wells (positive = injection)")
lines!(ax, t_jutul, wrat_cold, label = "COLD well", color = :blue)
lines!(ax, t_jutul, wrat_warm, label = "WARM well", color = :orange)
axislegend(position = :ct)
fig
# ## Plot the final temperature in the reservoir
# We see the final temperature distribution in the reservoir. The regions near
# both wells are warmer than the rest of the reservoir, with the hot well being
# the warmest.
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1], zreversed = true, azimuth = 4.55, elevation = 0.2)
for (wname, w) in get_model_wells(case)
    plot_well!(ax, G, w)
end
Jutul.plot_mesh_edges!(ax, G, alpha = 0.1)
temp = states[end][:Temperature] .- 273.15
plt = plot_cell_data!(ax, G, temp,
    colormap = :thermal,
    cells = findall(x -> x > 18, temp),
    transparency = true,
    shading = NoShading
)
Colorbar(fig[2, 1], plt, label = "Temperature (°C)", vertical = false)
fig
# ## Plot the temperature near the warm well and compare to E300
# We compare the temperature in cells close to the warm well in JutulDarcy and
# the same case simulated in E300, demonstrating excellent agreement.
warm1 = map(x -> x[:Temperature][cell_warm1] - 273.15, states)
warm2 = map(x -> x[:Temperature][cell_warm2] - 273.15, states)

t_e300 = cdata[:, 1]
warm1_e300 = cdata[:, 2]
warm2_e300 = cdata[:, 3]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time elapsed (days)", ylabel = "Temperature", title = "Temperature in cells close to warm well")
lines!(ax, t_jutul, warm1, label = "JutulDarcy, cell (105,75,6)", color = :orange)
lines!(ax, t_e300, warm1_e300, label = "E300, cell (105,75,6)", linestyle = :dash, linewidth = 3, color = :orange)

lines!(ax, t_jutul, warm2, label = "JutulDarcy, cell (106,75,6)", color = :blue)
lines!(ax, t_e300, warm2_e300, label = "E300, cell (106,75,6)", linestyle = :dash, linewidth = 3, color = :blue)
axislegend(position = :cb)
fig
# ## Plot the temperature near the cold well and compare
# We note a similar match between the solvers near the cold well.
cold1 = map(x -> x[:Temperature][cell_cold1] - 273.15, states)
cold2 = map(x -> x[:Temperature][cell_cold2] - 273.15, states)

t_e300 = cdata[:, 1]
cold1_e300 = cdata[:, 4]
cold2_e300 = cdata[:, 5]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time elapsed (days)", ylabel = "Temperature", title = "Temperature in cells close to warm well")
lines!(ax, t_jutul, cold1, label = "JutulDarcy, cell (50,75,6)", color = :orange)
lines!(ax, t_e300, cold1_e300, label = "E300, cell (50,75,6)", linestyle = :dash, linewidth = 3, color = :orange)
lines!(ax, t_jutul, cold2, label = "JutulDarcy, cell (51,75,6)", color = :blue)
lines!(ax, t_e300, cold2_e300, label = "E300, cell (51,75,6)", linestyle = :dash, linewidth = 3, color = :blue)
axislegend(position = :ct)
fig
# ## Plot the well temperatures
# Finally, we compare the reported temperatures in the wells between JutulDarcy
# and E300. These values are a mix of prescribed conditions (during injection)
# and the solution values (during production).
t_well_e300 = wdata[:, 1]
warm_e300 = wdata[:, 2]
warm_jutul = ws[:WARM][:temperature] .- 273.15
cold_e300 = wdata[:, 3]
cold_jutul = ws[:COLD][:temperature] .- 273.15

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time elapsed (days)", ylabel = "Temperature (°C)", title = "Well temperatures")
lines!(ax, t_jutul, warm_jutul, label = "JutulDarcy, WARM", color = :orange)
lines!(ax, t_well_e300, warm_e300, label = "E300, WARM", linestyle = :dash, linewidth = 3, color = :orange)
lines!(ax, t_jutul, cold_jutul, label = "JutulDarcy, COLD", color = :blue)
lines!(ax, t_well_e300, cold_e300, label = "E300, COLD", linestyle = :dash, linewidth = 3, color = :blue)
axislegend(position = :cb)
fig
# ## Plot the total energy in the reservoir
# We plot the total energy in the reservoir relative to the initial condition by
# summing up the thermal energy in all cells at the initial state and using this
# as the baseline. The total energy in the reservoir matches the injected and
# produced energy, as there are no open boundary conditions.
E0 = sum(states[1][:TotalThermalEnergy])
energy = map(x -> sum(x[:TotalThermalEnergy]) - E0, states)
lines(t_jutul, energy, axis = (xlabel = "Time elapsed (days)", ylabel = "Total energy (J)", title = "Total energy in reservoir relative to initial condition"))
# ## Plot the reservoir in the interactive viewer
# If you are running this example yourself, you can launch an interactive viewer
# and explore the evolution of the model.
plot_reservoir(case, states, key = :Pressure, step = 100)
