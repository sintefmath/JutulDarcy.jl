# # Gradient-based optimization of net present value (NPV)
# One usage of a differentiable simulator is control optimization. A counterpart
# to history matching, control optimization is the process of finding the
# optimal controls for a given objective and simulation model. In this example,
# we will optimize the injection rates for a simplified version of the EGG model
# to maximize the net present value (NPV) of the reservoir.

# ## Setting up a coarse model
# We start off by loading the Egg model and coarsening it to a 20x20x3 grid. We
# limit the optimization to the first 50 timesteps to speed up the optimization
# process. If you want to run the optimization for all timesteps on the fine
# model directly, you can remove the slicing of the case and replace the coarse
# case with the fine case.
using Jutul, JutulDarcy, GLMakie, GeoEnergyIO, HYPRE, LBFGSB
data_dir = GeoEnergyIO.test_input_file_path("EGG")
data_pth = joinpath(data_dir, "EGG.DATA")
fine_case = setup_case_from_data_file(data_pth)
fine_case = fine_case[1:50];
coarse_case = coarsen_reservoir_case(fine_case, (20, 20, 3), method = :ijk);
# ## Set up the rate optimization
# We use an utility to set up the rate optimization problem. The utility sets up
# the objective function, constraints, and initial guess for the optimization.
# The optimization is set up to maximize the NPV of the reservoir. The
# contribution to the NPV for a given time ``t_i`` given in years is defined as:
# ``\text{NPV}_i = \Delta t_i C_i(q_o, q_w)(1 + r)^{-t_i}``
#
# Here, ``r`` is the discount rate and ``C_i(q_o, q_w)`` is the cash flow for
# the current rates. The cash flow is defined as the price of producing each
# phase (oil and water, with oil having a positive price and water a negative
# price) minus the cost of injecting water.
#
# We set the prices and costs (per barrel) as well as a discount rate of 5% per
# year. The base rate will be used for all injectors initially, which matches
# the base case for the Egg model. The utility function takes in a `steps`
# argument that can be used to set how often the rates are allowed to change.
#
# The values used here are arbitrary for the purposes of the example. You are
# encouraged to play around with the values to see how the outcome changes in
# the optimization. Generally higher discount rates prioritze immediate income,
# while lower or zero discount rates prioritize maximizing the oil production
# over the entire simulation period.
ctrl = coarse_case.forces[1][:Facility].control
base_rate = ctrl[:INJECT1].target.value
function optimize_rates(steps; use_box_bfgs = true)
    setup = JutulDarcy.setup_rate_optimization_objective(coarse_case, base_rate,
        max_rate_factor = 10,
        oil_price = 100.0,
        water_price = -10.0,
        water_cost = 5.0,
        discount_rate = 0.05,
        maximize = use_box_bfgs,
        sim_arg = (
            rtol = 1e-5,
            tol_cnv = 1e-5
        ),
        steps = steps
    )
    if use_box_bfgs
        obj_best, x_best, hist = Jutul.unit_box_bfgs(setup.x0, setup.obj,
            maximize = true,
            lin_eq = setup.lin_eq
        )
        H = hist.val
    else
        lower = zeros(length(setup.x0))
        upper = ones(length(setup.x0))
        results, x_best = lbfgsb(setup.F!, setup.dF!, setup.x0,
            lb=lower,
            ub=upper,
            iprint = 1,
            factr = 1e12,
            maxfun = 20,
            maxiter = 20,
            m = 20
        )
        H = results
    end
    return (setup.case, H, x_best)
end
# ## Optimize the rates
# We optimize the rates for two different strategies. The first strategy is to
# optimize constant rates for the entire period. The second strategy is to
# optimize the rates for each report step.

# ### Optimize with a single set of rates
case1, hist1, x1 = optimize_rates(:first);
# ### Plot rate allocation for the constant case
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Injector number", ylabel = "Rate fraction (of max injection rate)")
barplot!(ax, x1)
ax.xticks = eachindex(x1)
fig
# ### Optimize with varying rates per time-step
case2, hist2, x2 = optimize_rates(:each);
# ### Plot rate allocation per well
allocs = reshape(x2, length(x1), :)
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Report step", ylabel = "Rate fraction (of max injection rate)")
for i in axes(allocs, 1)
    lines!(ax, allocs[i, :], label = "Injector #$i")
end
axislegend()
fig
# ## Plot the evolution of the NPV
# We plot the evolution of the NPV for the two strategies to compare the
# results. Note that the optimization produces a higher NPV for the varying
# rates. This is expected, as the optimization can adjust the rates to the
# changes in the mobility field during the progress of the simulation.
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "LBFGS iteration", ylabel = "Net present value (million USD)")
scatter!(ax, 1:length(hist1), hist1./1e6, label = "Constant rates")
scatter!(ax, 1:length(hist2), hist2./1e6, marker = :x, label = "Varying rates")
axislegend(position = :rb)
fig
# ## Simulate the results
# We finally simulate all the cases to compare the results.
ws0, states0 = simulate_reservoir(coarse_case, info_level = -1)
ws1, states1 = simulate_reservoir(case1, info_level = -1)
ws2, states2 = simulate_reservoir(case2, info_level = -1)
# ### Compute measurables and compare
# NPV favors early production, as later income/costs have less relative value.
# We compute the field measurables to be able to compare the oil and water
# production.
f0 = reservoir_measurables(coarse_case.model, ws0)
f1 = reservoir_measurables(case1.model, ws1)
f2 = reservoir_measurables(case2.model, ws2)
# ### Plot the cumulative field oil production
# The optimized cases produce more oil than the base case.
bbl = si_unit(:stb)
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time / days", ylabel = "Field oil production (accumulated, barrels)", title = "Base case")
t = ws0.time./si_unit(:day)
lines!(ax, t, f0[:fopt].values./bbl, label = "Base case")
lines!(ax, t, f1[:fopt].values./bbl, label = "Constant rates")
lines!(ax, t, f2[:fopt].values./bbl, label = "Varying rates")
axislegend(position = :rb)
fig
# ### Plot the cumulative field water production
# The optimized cases produce less water than the base case.
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time / days", ylabel = "Field water production (accumulated, barrels)", title = "Base case")
t = ws0.time./si_unit(:day)
lines!(ax, t, f0[:fwpt].values./bbl, label = "Base case")
lines!(ax, t, f1[:fwpt].values./bbl, label = "Constant rates")
lines!(ax, t, f2[:fwpt].values./bbl, label = "Varying rates")
axislegend(position = :rb)
fig
# ### Plot the differences in water saturation
# The optimized cases have a different water saturation distribution compared to the base case.
reservoir = reservoir_domain(coarse_case.model)
g = physical_representation(reservoir)
function plot_diff!(ax, s)
    plt = plot_cell_data!(ax, g, s0 - s1, colorrange = (-0.5, 0.5), colormap = :balance)
    for (w, wd) in get_model_wells(coarse_case.model)
        plot_well!(ax, g, wd, top_factor = 0.8, fontsize = 10)
    end
    ax.elevation[] = 1.03
    ax.azimuth[] = 3.75
    return plt
end

s0 = states0[end][:Saturations][1, :]
s1 = states1[end][:Saturations][1, :]
s2 = states2[end][:Saturations][1, :]
fig = Figure(size = (1600, 800))
ax = Axis3(fig[1, 1], zreversed = true, title = "Water saturation difference, constant rates")
plt = plot_diff!(ax, s1)
ax = Axis3(fig[1, 2], zreversed = true, title = "Water saturation difference, varying rates")
plot_diff!(ax, s2)
Colorbar(fig[2, :], plt, vertical = false)
fig