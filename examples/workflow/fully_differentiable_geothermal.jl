# # A fully differentiable geothermal doublet: History matching and control optimization
# We are going to set up a conceptual geothermal doublet model in 2D and perform
# gradient based history matching. This example serves two main purposes:
# 1. It demonstrates the conceptual workflow for setting up a geothermal model
#    from scratch with a fairly straightforward mesh setup.
# 2. It shows how to set up a gradient based history matching workflow with the
#    generic optimization interface that allows for optimizing any input
#    parameter used in the setup of a model.
# ## Load packages and define units
using Jutul, JutulDarcy, HYPRE, GeoEnergyIO, GLMakie
meter, kilogram, bar, year, liter, second, darcy, day = si_units(:meter, :kilogram, :bar, :year, :liter, :second, :darcy, :day)

# ## Set up the reservoir mesh
# The model is a typical geothermal case where there is a layer of high
# permeability in the middle, confined between two low-permeable layers. For a
# geothermal model, the low permeable layers are important, as they store
# significant amounts of heat that can be conducted to the high permeable layer
# during production.
#
# We set up the mesh so that the high permeable layer where most of the
# advective transport occurs has a higher lateral resolution than the low
# permeable layers. The model is also essentally 2D as there is only one cell
# thickness in the y direction - a choice that is made to make the example fast
# to run, especially during the later optimization stages where many simulations
# must be run to achieve convergence.
nx = 50
ntop = 5
nmiddle = 10
nbottom = 5
nz = ntop + nmiddle + nbottom
# ### Set up layer thicknesses and vertical cell thicknesses
top_layer_thickness = 300.0*meter
middle_layer_thickness = 200.0*meter
bottom_layer_thickness = 300.0*meter
dz = Float64[]
for i in 1:ntop
    push!(dz, top_layer_thickness/ntop)
end
for i in 1:nmiddle
    push!(dz, middle_layer_thickness/nmiddle)
end
for i in 1:ntop
    push!(dz, bottom_layer_thickness/nbottom)
end

cmesh = CartesianMesh((nx, 1, nz), (2000.0, 50.0, dz))
rmesh = UnstructuredMesh(cmesh, z_is_depth = true)
# ### Define regions based on our selected depths
# We tag each cell with a region number based on its depth. The top layer is
# region 1, the middle layer is region 2, and the bottom layer is region 3.
geo = tpfv_geometry(rmesh)
depths = geo.cell_centroids[3, :]
regions = Int[]
for (i, d_i) in enumerate(depths)
    if d_i <= top_layer_thickness
        r = 1
    elseif d_i <= top_layer_thickness + middle_layer_thickness
        r = 2
    else
        r = 3
    end
    push!(regions, r)
end
# ### Plot the mesh and regions
fig, ax, plt = plot_cell_data(rmesh, regions,
    alpha = 0.5,
    outer = true,
    transparency = true,
    colormap = Categorical(:heat)
)
ax.elevation[] = 0.0
ax.azimuth[] = π/2
plot_mesh_edges!(ax, rmesh)
fig
# ## Define functions for setting up the simulation
# We will define a function that takes in a Dict with different values and sets
# up the simulation. The key idea is that we can then optimize the values in the
# Dict to perform optimization. As we can define any such Dict to set up the
# model, this interface is very flexible and can be used for both control
# optimization and history matching with respect to almost any parameter of the
# model. The disadvantage is that the setup function will be called many times,
# which can be a substantial cost compared to the more structured optimization
# interface that only allows for optimization of the numerical parameters (e.g.
# for the CGNet example).

# ### Define the time schedule
# We set up a time schedule for the simulation. The total simulation time is 30
# years, and we report the results every 120 days. We also define ten different
# intervals in this 30 year period, which are the period where we will allow the
# rates and temperatures to vary during the last part of the optimization tutorial.
total_time = 30.0*year
report_step_length = 120.0*day
dt = fill(report_step_length, Int(ceil(total_time/report_step_length)))
num_intervals = 10
interval_interval = total_time/num_intervals
interval_for_step = map(t -> min(Int(ceil(t/interval_interval)), num_intervals), cumsum(dt))
# ### Define the wells
# We set up two wells, one injector and one producer. The injector is located at
# the left side of the model, and the producer is located at the right side. We
# use multisegment wells.
base_rate = 15*liter/second
base_temp = 15.0

domain = reservoir_domain(rmesh)
inj_well = setup_vertical_well(domain, 5, 1,
    heel = ntop+1,
    toe = ntop+nmiddle,
    name = :Injector,
    simple_well = false
)
prod_well = setup_vertical_well(domain, nx - 5, 1,
    heel = ntop+1,
    toe = ntop+nmiddle,
    name = :Producer,
    simple_well = false
)

model_base, = setup_reservoir_model(
    domain, :geothermal,
    wells = [inj_well, prod_well],
);
# ### Set up a helper to define the forces for a given rate and temperature
function setup_doublet_forces(model, inj_temp, inj_rate)
    T_Kelvin = convert_to_si(inj_temp, :Celsius)
    rate_target = TotalRateTarget(inj_rate)
    ctrl_inj  = InjectorControl(rate_target, [1.0],
        density = 1000.0, temperature = T_Kelvin)

    bhp_target = BottomHolePressureTarget(50*bar)
    ctrl_prod = ProducerControl(bhp_target)

    control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)
    return setup_reservoir_forces(model, control = control)
end
# ### Define the main setup function
# This function sets up the model based on the parameters provided in the Dict.
# It takes in two arguments: The required parameters in a Dict and an optional
# step_info argument that can be used to set up the model for a specific time
# step. The function returns a JutulCase object that can be used to simulate the
# reservoir. Here, we ignore the `step_info` argument and set up the entire
# schedule every time. Jutul will then automatically use the correct force based
# on the time step in the simulation.

function setup_doublet_case(prm, step_info = missing)
    model = deepcopy(model_base)
    rdomain = reservoir_domain(model)
    rdomain[:permeability] = prm["layer_perm"][regions]
    rdomain[:porosity] = prm["layer_porosities"][regions]
    rdomain[:rock_heat_capacity] = prm["layer_heat_capacity"][regions]

    T0 = convert_to_si(70, :Celsius)
    thermal_gradient = 20.0/1000.0*meter
    eql = EquilibriumRegion(model, 50*bar, 0.0, temperature_vs_depth = z -> T0 + z*thermal_gradient)
    state0 = setup_reservoir_state(model, eql)

    forces_per_interval = map((T, rate) -> setup_doublet_forces(model, T, rate),
        prm["injection_temperature_C"], prm["injection_rate"])

    forces = forces_per_interval[interval_for_step]

    return JutulCase(model, dt, forces, state0 = state0)
end

# ## Perform a history match
# We first set up a truth case that we will use to generate the data for the
# history match. We define high perm and porosity in the middle layer, and low
# perm and porosity in the top and bottom layers before simulating the model.
prm_truth = Dict(
    "injection_rate" => fill(base_rate, num_intervals),
    "injection_temperature_C" => fill(base_temp, num_intervals),
    "layer_porosities" => [0.1, 0.3, 0.1],
    "layer_perm" => [0.01, 0.8, 0.02].*darcy,
    "layer_heat_capacity" => [500.0, 600.0, 450.0], # Watt / m K
)
case_truth = setup_doublet_case(prm_truth)
ws, states = simulate_reservoir(case_truth)
# ### Define a mismatch objective function
# The mismatch objective function is defined as the sum of squares difference
# between the simulated values and the reference values observed in the wells.
# Note that we only make use of the well data:
# - The temperature at the producer well
# - The mass rate at the producer well (since it is controlled on BHP)
# - The BHP at the injector well (since it is controlled on rate)
#
# We use the `get_1d_interpolator` function to create interpolators for the
# reference values, since we cannot assume that the simulator will use exactly
# the same time-steps as the reference values.
prod_rate = ws.wells[:Producer][:mass_rate]
prod_temp = ws.wells[:Producer][:temperature]
inj_bhp = ws.wells[:Injector][:bhp]

prod_temp_by_time = get_1d_interpolator(ws.time, prod_temp)
prod_rate_by_time = get_1d_interpolator(ws.time, prod_rate)
inj_pressure_by_time = get_1d_interpolator(ws.time, inj_bhp)

import JutulDarcy: compute_well_qoi
function mismatch_objective(m, s, dt, step_info, forces)
    current_time = step_info[:time] + dt
    ## Current values
    T_at_prod = compute_well_qoi(m, s, forces, :Producer, :temperature)
    mass_rate = compute_well_qoi(m, s, forces, :Producer, :mass_rate)
    bhp = compute_well_qoi(m, s, forces, :Injector, :bhp)
    ## Reference values
    T_at_prod_ref = prod_temp_by_time(current_time)
    mass_rate_ref = prod_rate_by_time(current_time)
    bhp_ref = inj_pressure_by_time(current_time)
    ## Define mismatch by scaling each term
    T_mismatch = (T_at_prod_ref - T_at_prod)*100.0
    rate_mismatch = (mass_rate_ref - mass_rate)
    bhp_mismatch = (bhp - bhp_ref)/(50*bar)
    return dt * sqrt(T_mismatch^2 + rate_mismatch^2 + bhp_mismatch^2) / total_time
end
# ### Pick an initial guess
# We set up an initial guess for the parameters that we will optimize. We assume
# the injection rate and temperature to be known and we set the porosities and
# permeabilities to uniform values. The heat capacity is given a bit of
# layering, but still with completely wrong values.
prm_guess = Dict(
    "injection_rate" => fill(base_rate, num_intervals),
    "injection_temperature_C" => fill(base_temp, num_intervals),
    "layer_porosities" => [0.2, 0.2, 0.2],
    "layer_perm" => [0.2, 0.2, 0.2].*darcy,
    "layer_heat_capacity" => [300.0, 200.0, 400.0]
)
case_guess = setup_doublet_case(prm_guess)
ws_guess, states_guess = simulate_reservoir(case_guess)
# ### Set up the optimization
# We define a dictionary optimization problem that will optimize the parameters
# in the `prm_guess` dictionary. We start by setting up the object itself, which
# takes in the initial guess `Dict` and the corresponding setup function.
opt = JutulDarcy.setup_reservoir_dict_optimization(prm_guess, setup_doublet_case)
# ### Define active parameters and their limits
# Note that while the parameters get listed, they are all marked as inactive. We
# need to explicitly make them free/active and specify a range for each
# parameter before we can optimize them. We use wide absolute limits for each entry.
free_optimization_parameter!(opt, "layer_perm", abs_max = 1.5*darcy, abs_min = 0.01*darcy)
free_optimization_parameter!(opt, "layer_heat_capacity", abs_max = 1000.0, abs_min = 100.0)
free_optimization_parameter!(opt, "layer_porosities", abs_max = 0.5, abs_min = 0.05)
# ### Call the optimizer
# Now that we have freed a few parameters, we can call the optimizer with the
# objective function. The defaults for the optimizer are fairly reasonable, so
# we do not tweak the convergence criteria or the maximum number of iterations.
# Note that by default the optimizer uses LBFGS, but it is also possible to pass
# other optimizers as a function callable. The default for the optimizer is to
# minimize the objective function, which is the case for a history match. By
# passing for example `lbfgs_num = 1, max_it = 50` it is possible to obtain a
# better match, but this is not necessary for the purpose of this example.
prm_opt = JutulDarcy.optimize_reservoir(opt, mismatch_objective);
# ### Print the optimization overview
# If we display the optimization overview, we can see that there are now
# additional columns indicating the optimized values. Note that while the
# permeability and porosities are well matched, the heat capacity of the low
# permeable layers are not very accurate. There is likely not enough data in the
# production profiles to constrain the heat capacity of the low permeable
# layers, as there is limited heat siphoned from these layers in the truth case.
opt
# ### Simulate the optimized case
case_opt = setup_doublet_case(prm_opt)
ws_opt, states_opt = simulate_reservoir(case_opt)
# ### Plot the results
# We plot the results for the truth case, the initial guess, and the optimized
# case. The temperature is plotted in Celsius and we use the same color scale
# for all steps.
step = 80
cmap = :heat
fig = Figure(size = (1200, 400))
ax = Axis3(fig[1, 1], title = "Truth")
plot_cell_data!(ax, rmesh, states[step][:Temperature] .- 273.15, colorrange = (10.0, 100.0), colormap = cmap)
ax.elevation[] = 0.0
ax.azimuth[] = -π/2
hidedecorations!(ax)

ax = Axis3(fig[1, 2], title = "Initial guess")
plot_cell_data!(ax, rmesh, states_guess[step][:Temperature] .- 273.15, colorrange = (10.0, 100.0), colormap = cmap)
ax.elevation[] = 0.0
ax.azimuth[] = -π/2
hidedecorations!(ax)

ax = Axis3(fig[1, 3], title = "Optimized")
plt = plot_cell_data!(ax, rmesh, states_opt[step][:Temperature] .- 273.15, colorrange = (10.0, 100.0), colormap = cmap)
ax.elevation[] = 0.0
ax.azimuth[] = -π/2
hidedecorations!(ax)
Colorbar(fig[2, 1:3], plt, vertical = false)
fig
# ## Set up control optimization
# We can also use the same setup to perform control optimization, where we now
# can take advantage of the per-interval selection of rates and temperatures.
# Admittely, this problems is fairly simple, so the optimization is more
# conceptual than realistic: We define a new objective function that uses a
# fixed cost for the injected water (per degree times rate) and a similar value
# of produced heat. To make the optimization problem non-trivial, the cost of
# additional water (or higher temperature water) is significantly higher than
# the value of produced water with the same temperature.
temperature_injection_cost = 15.0
temperature_production_value = 8.0

function optimization_objective(m, s, dt, step_info, forces)
    current_time = step_info[:time] + dt
    T_at_prod = convert_from_si(compute_well_qoi(m, s, forces, :Producer, :temperature), :Celsius)
    T_at_inj = convert_from_si(forces[:Facility].control[:Injector].temperature, :Celsius)

    mass_rate_injector = compute_well_qoi(m, s, forces, :Injector, :mass_rate)
    mass_rate_producer = compute_well_qoi(m, s, forces, :Producer, :mass_rate)

    cost_inj = abs(mass_rate_injector) * T_at_inj * temperature_injection_cost
    value_prod = abs(mass_rate_producer) * T_at_prod * temperature_production_value
    return dt * (value_prod - cost_inj) / total_time
end

opt_ctrl = JutulDarcy.setup_reservoir_dict_optimization(prm_truth, setup_doublet_case)
# ### Set optimization to use injection rate and temperature
# Note that as these are represented as per-interval values, we could also have passed vectors of equal length as the
# number of intervals for more fine-grained control over the limits.
free_optimization_parameter!(opt_ctrl, "injection_temperature_C", abs_max = 80.0, abs_min = 10.0)
free_optimization_parameter!(opt_ctrl, "injection_rate", abs_min = 1.0*liter/second, abs_max = 30.0*liter/second)
# ### Call the optimizer
prm_opt_ctrl = JutulDarcy.optimize_reservoir(opt_ctrl, optimization_objective, maximize = true);
opt_ctrl
# ### Plot the optimized injection rates and temperatures
# The optimized injection rates and temperatures are plotted for each interval.
# The base case is shown in blue, while the optimized case is shown in orange.
# Note that the optimized case has reduced the injection temperature to the
# lower limit for all steps, and instead increase the injection rate
# significantly. The injection rate has a decrease part-way during the
# simulation, which increases the residence time of the injected water, allowing
# additional heat to be siphoned from the low permeable layers.
fig = Figure(size = (1200, 400))
ax = Axis(fig[1, 1], title = "Optimized injection temperature", ylabel = "Injection temperature (°C)", xlabel = "Interval")
scatter!(ax, prm_truth["injection_temperature_C"], label = "Base case")
scatter!(ax, prm_opt_ctrl["injection_temperature_C"], label = "Optimized case")
axislegend(position = :rc)

ax = Axis(fig[1, 2], title = "Optimized injection rate", ylabel = "Liter/second", xlabel = "Interval")
scatter!(ax, prm_truth["injection_rate"]./(liter/second), label = "Base case")
scatter!(ax, prm_opt_ctrl["injection_rate"]./(liter/second), label = "Optimized case")
axislegend(position = :rc)
fig
# ### Simulate the optimized case
case_opt_ctrl = setup_doublet_case(prm_opt_ctrl)
ws_opt_ctrl, states_opt_ctrl = simulate_reservoir(case_opt_ctrl)
# ### Plot the distribution of temperature with and without optimization
step = 80
cmap = :heat
fig = Figure(size = (1000, 400))
ax = Axis3(fig[1, 1], title = "Base case")
plot_cell_data!(ax, rmesh, states[step][:Temperature] .- 273.15, colorrange = (10.0, 100.0), colormap = cmap)
ax.elevation[] = 0.0
ax.azimuth[] = -π/2
hidedecorations!(ax)

ax = Axis3(fig[1, 2], title = "Optimized")
plt = plot_cell_data!(ax, rmesh, states_opt_ctrl[step][:Temperature] .- 273.15, colorrange = (10.0, 100.0), colormap = cmap)
ax.elevation[] = 0.0
ax.azimuth[] = -π/2
hidedecorations!(ax)
Colorbar(fig[2, 1:2], plt, vertical = false)
fig
# ### Plot the total thermal energy in the reservoir
# The total thermal energy in the reservoir is computed as the sum of the
# thermal energy in each cell, which is the result of the rock heat capacity,
# porosity, fluid heat capacity and the temperature in each cell. The optimized
# strategy significantly decreases the remaining thermal energy in the
# reservoir, while still producing less cost than the base case according to our
# objective. The 2D nature of this problem makes it easy to recover a large
# amount energy, as the majority of the cells are swept by the cold front.
total_energy = map(s -> sum(s[:TotalThermalEnergy]), states)
total_energy_opt = map(s -> sum(s[:TotalThermalEnergy]), states_opt_ctrl)

fig = Figure(size = (1200, 400))
ax = Axis(fig[1, 1], title = "Total thermal energy", ylabel = "Total remaining energy (megajoules)", xlabel = "Time (days)")
t = ws.time ./ si_unit(:day)
lines!(ax, t, total_energy./1e6, label = "Base case")
lines!(ax, t, total_energy_opt./1e6, label = "Optimized case")
axislegend(position = :rc)
fig