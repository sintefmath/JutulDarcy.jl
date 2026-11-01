# # Geothermal well groups
# <tags: Geothermal, Wells, Introduction> Three injectors maintain pressure
# while four producers supply hot water. Group targets share a surface rate that
# they collectively aim to achieve, with individual wells subject to their own
# limits (e.g. on bottom-hole-pressure).
using Jutul, JutulDarcy

bar, day, year, darcy = si_units(:bar, :day, :year, :darcy)

# ## Reservoir and wells
# A smooth permeability variation gives the wells different injectivities.
mesh = CartesianMesh((15, 11, 2), (1200.0, 800.0, 60.0))
geometry = tpfv_geometry(mesh)
permeability = [0.25darcy*exp(0.6*sin(2pi*x/1200)*cos(2pi*y/800))
    for (x, y, z) in eachcol(geometry.cell_centroids)]
domain = reservoir_domain(mesh; permeability, porosity = 0.18)

injectors = [:I1, :I2, :I3]
producers = [:P1, :P2, :P3, :P4]
wells = [setup_vertical_well(domain, 3, j; name, simple_well = true)
    for (name, j) in zip(injectors, (3, 6, 9))]
append!(wells, [setup_vertical_well(domain, 13, j; name, simple_well = true)
    for (name, j) in zip(producers, (2, 5, 8, 10))])

# ## Group allocation
# Use equal weights for injection and the default potential weights for production.
# Every well in a group receives the same target for the entire group.
model = setup_reservoir_model(domain, :geothermal;
    wells,
    well_groups = Dict(:Injection => injectors, :Production => producers),
    group_use_potentials = Dict(:Injection => false)
)
state0 = setup_reservoir_state(model;
    Pressure = 120bar, Temperature = convert_to_si(140.0, :Celsius))

# Both groups request 1200 surface cubic metres per day of water.
# Injection is positive and production is negative.
rate = 1200.0/day
injection = InjectorControl(GroupSurfaceWaterRateTarget(rate), [1.0];
    density = 1000.0, temperature = convert_to_si(40.0, :Celsius))
production = ProducerControl(GroupSurfaceWaterRateTarget(-rate))
control = Dict{Symbol, Any}(w => injection for w in injectors)
merge!(control, Dict(w => production for w in producers))

# ## Individual limits
# Pressure limits protect all wells. I1 and P1 also have tighter rate limits;
# their unused shares are redistributed among the other group members.
limits = Dict{Symbol, Any}(w => (bhp = 170bar,) for w in injectors)
merge!(limits, Dict(w => (bhp = 70bar,) for w in producers))
limits[:I1] = (bhp = 170bar, wrat = 250.0/day)
limits[:P1] = (bhp = 70bar, wrat = -180.0/day)
forces = setup_reservoir_forces(model; control, limits)

# ## Simulate and inspect
# Run ten years, with shorter initial steps to establish the well controls.
dt = [5day, 10day, 15day, fill((10year - 30day)/40, 40)...]
result = simulate_reservoir(state0, model, dt; forces)

using GLMakie
plot_reservoir(model, result.states)
plot_well_results(result.wells)

# Group totals stay at their targets while the remaining wells have capacity.
# Production eventually falls below its target as producers reach the BHP limit.
# Compare these curves with individual rates and BHPs in the well plot.
time = result.time/year
injected = sum(result.wells[w][:wrat] for w in injectors)*day
produced = -sum(result.wells[w][:wrat] for w in producers)*day
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Time [years]", ylabel = "Water rate [m³/day]")
lines!(ax, time, injected; label = "Injection group")
lines!(ax, time, produced; label = "Production group")
hlines!(ax, [rate*day]; color = :black, linestyle = :dash, label = "Group target")
axislegend(ax)
fig
