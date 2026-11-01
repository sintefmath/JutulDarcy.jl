# # Geothermal well groups
# <tags: Geothermal, Wells, Introduction> Three injectors maintain pressure
# while four producers supply hot water. Group targets share a surface rate that
# they collectively aim to achieve, with individual wells subject to their own
# limits (e.g. on bottom-hole-pressure).
using Jutul, JutulDarcy, GLMakie
bar, day, year, darcy = si_units(:bar, :day, :year, :darcy)

# ## Reservoir and wells
# A smooth permeability variation gives the wells a bit different rates as the
# pressure limit and consequently different allocations of the group rate.
mesh = CartesianMesh((15, 11, 2), (1200.0, 800.0, 60.0))
geometry = tpfv_geometry(mesh)
permeability = zeros(number_of_cells(mesh))
for cell in 1:number_of_cells(mesh)
    x = geometry.cell_centroids[1, cell]
    y = geometry.cell_centroids[2, cell]
    permeability[cell] = 0.25darcy*exp(0.6*sin(2pi*x/1200)*cos(2pi*y/800))
end
domain = reservoir_domain(mesh; permeability, porosity = 0.18)

injectors = [:I1, :I2, :I3]
producers = [:P1, :P2, :P3, :P4]
injector_rows = [3, 6, 9]
producer_rows = [2, 5, 8, 10]
wells = []
for i in eachindex(injectors)
    well = setup_vertical_well(domain, 3, injector_rows[i];
        name = injectors[i], simple_well = true)
    push!(wells, well)
end
for i in eachindex(producers)
    well = setup_vertical_well(domain, 13, producer_rows[i];
        name = producers[i], simple_well = true)
    push!(wells, well)
end

# ## Group allocation
# Group allocation uses potential weights by default where each well receives a
# share based on its maximal potential for injection or production at the
# corresponding BHP limit. Set a group to false for equal weights. Every well
# in a group receives the same target for the entire group.
#
# We also disable shutting of wells, since this can lead to multiple solutions
# for group control cases.
model = setup_reservoir_model(domain, :geothermal;
    wells,
    well_groups = Dict(:Injection => injectors, :Production => producers),
    group_use_potentials = Dict(:Injection => false),
    can_shut_wells = false
)

# Initialize hydrostatic pressure from 120 bar at the reservoir's mid-depth.
# One equilibrium region covers the reservoir, initially at 140 degrees Celsius.
equilibrium = EquilibriumRegion(model, 120bar, 30.0;
    temperature = convert_to_si(140.0, :Celsius))
state0 = setup_reservoir_state(model, equilibrium)

## Inspect permeability and well locations
plot_reservoir(model; key = :permeability, fancy = false)

# ## Set up well constraints
# Both groups request 1200 surface cubic metres per day of water. Injection is
# positive and production is negative. We set the same group rate for all wells
# in a group.
rate = 1200.0/day
injection = InjectorControl(GroupSurfaceWaterRateTarget(rate), [1.0];
    density = 1000.0, temperature = convert_to_si(40.0, :Celsius))
production = ProducerControl(GroupSurfaceWaterRateTarget(-rate))
control = Dict{Symbol, Any}()
for well in injectors
    control[well] = injection
end
for well in producers
    control[well] = production
end

# ## Individual limits
# Pressure limits protect all wells. I1 and P1 also have tighter rate limits;
# their unused shares are redistributed among the other group members when they
# are operating on the limits.
limits = Dict{Symbol, Any}()
for well in injectors
    limits[well] = (bhp = 170bar,)
end
for well in producers
    limits[well] = (bhp = 70bar,)
end
limits[:I1] = (bhp = 170bar, wrat = 250.0/day)
limits[:P1] = (bhp = 70bar, wrat = -180.0/day)
forces = setup_reservoir_forces(model; control, limits)

# ## Simulate and inspect
# Run ten years, with shorter initial steps to establish the well controls.
dt = [5day, 10day, 15day]
for step in 1:40
    push!(dt, (10year - 30day)/40)
end
dt = fill(31day, 12*10)
result = simulate_reservoir(state0, model, dt; forces = forces, output_substates = true)

# Show the final temperature distribution and the individual well histories.
plot_reservoir(model, result.states; key = :Temperature, step = length(result.states), fancy = false)
## Plot group summaries
# Group totals stay at their targets while the remaining wells have capacity.
# Production eventually falls below its target as producers reach the BHP limit.
# Compare these curves with individual rates and BHPs in the well plot.
# Group rates and cumulative volumes are included in the simulation summary.
plot_summary(result; plots = [
    "Injection:GWIR", "Production:GWPR",
    "Injection:GWIT", "Production:GWPT",
], cols = 2)
## Plot the individual producers to see the limits
plot_summary(result; plots = [
    "P1:WWPT", "P1:WBHP",
    "P2:WWPT", "P2:WBHP",
    "P3:WWPT", "P3:WBHP"
], cols = 3)
