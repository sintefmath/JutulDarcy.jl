# # Well control optimization for NPV (Egg model)
# <tags: Immiscible, InputFile, ModelReduction, Optimization, Advanced>
# This example is a JutulDarcy port of the MRST example
# [`optimizeEggEnsemble`](https://github.com/SINTEF-AppliedCompSci/MRST/blob/main/autodiff/optimization/examples/ensembles/optimizeEggEnsemble.m),
# restricted to a single realization. We optimize the well controls of the Egg
# model over a sequence of control periods to maximize net present value (NPV).
#
# Relative to [the rate optimization example](rate_optimization.md), which tunes
# a single injector rate per report step with a hand-built objective, this one
# uses the general [`setup_well_control_optimization`](@ref) helper: an explicit
# sequence of control periods, and an *extended control vector* that mixes well
# *targets* (injector rates, producer BHPs) and operating *limits* (an injector
# BHP cap) in the same optimization. The adjoint follows whichever control is
# active during the simulation, so an operating limit contributes to the
# gradient exactly on the steps where it binds and is inert otherwise. This
# mirrors MRST's generalized-reduced-gradient setup, which "computes gradients
# both for well target controls and well upper/lower simulation limits".

using Jutul, JutulDarcy, GLMakie, GeoEnergyIO, LinearAlgebra, Statistics, Printf

# ## Load the Egg model
# We load the first Egg realization from its `.DATA` file, keep the first five
# years, and coarsen to a 20x20x3 grid to keep the example fast. Replace the
# coarsening line with `egg_case = fine_case[1:n5]` to optimize on the full grid.
egg_dir = GeoEnergyIO.test_input_file_path("EGG")
fine_case = setup_case_from_data_file(joinpath(egg_dir, "EGG.DATA"))
n5 = count(<=(5*si_unit(:year)), cumsum(fine_case.dt))
egg_case = coarsen_reservoir_case(fine_case[1:n5], (20, 20, 3), method = :ijk)
nstep = length(egg_case.dt)

# ## The base schedule
# The Egg schedule has eight water injectors under total-rate control
# (79.5 m³/day, with a 420 bar BHP cap) and four producers under BHP control
# (395 bar). There is a single control that lasts the whole horizon.
injectors = [Symbol("INJECT$i") for i in 1:8]
producers = [Symbol("PROD$i")   for i in 1:4]

day  = si_unit(:day)
bar  = si_unit(:bar)
stb  = si_unit(:stb)

# ## Control periods
# Following MRST, we split the horizon into ten control periods of roughly equal
# length. Each well quantity gets an independent value in each period.
edges = round.(Int, range(0, nstep, length = 11))
control_periods = [(edges[i] + 1):edges[i + 1] for i in 1:10]

# ## The extended control vector
# `setup_well_control_optimization` takes a list of degree-of-freedom specs.
# A quantity that matches the well's base control target is optimized as the
# *target*; any other quantity is treated as an operating *limit* and written
# into `forces[:Facility].limits`. Bounds may be absolute (`abs_min`/`abs_max`)
# and/or relative to the base value (`rel_min`/`rel_max`).
#
# We use the bounds from the MRST example:
#
# | quantity | role | bounds |
# |---|---|---|
# | injector `rate` | target | `[1, 150]` m³/day |
# | injector `bhp` | limit (upper) | `[400, 420]` bar |
# | producer `bhp` | target | `[380, 400]` bar |
#
# A quantity is treated as the target when it matches the well's base control
# (injectors are rate-controlled, producers are BHP-controlled), so the injector
# `bhp` entry becomes an operating limit. A limit only binds on some steps of
# some evaluations: where it does not bind its gradient is zero and the optimizer
# leaves it alone. In this economic setting the injectors stay well below their
# BHP cap, so those entries stay put — but the adjoint does follow control
# switching that does happen (producers shut in, see below).
controls = NamedTuple[]
for w in injectors
    push!(controls, (well = w, quantity = :rate, abs_min = 1.0/day,  abs_max = 150.0/day))
    push!(controls, (well = w, quantity = :bhp,  abs_min = 400*bar,   abs_max = 420*bar))
end
for w in producers
    push!(controls, (well = w, quantity = :bhp,  abs_min = 380*bar,   abs_max = 400*bar))
end

# ## The NPV objective
# We reuse [`JutulDarcy.npv_objective`](@ref) with the economics from the MRST
# example: oil is worth 50 \$/stb, injected and produced water each cost
# 3 \$/stb. A 10 %/year discount rate rewards earlier production. `scale` brings
# the objective to order unity for the optimizer.
npv_arg = (
    oil_price     = 50.0,
    water_price   = -3.0,   # produced water: a cost
    water_cost    =  3.0,   # injected water: a cost
    oil_cost      =  0.0,
    discount_rate =  0.1,
    scale         =  1e8,
)

copt = setup_well_control_optimization(egg_case, control_periods, controls;
    objective = :npv, npv_arg = npv_arg)

# ## Verify the gradient
# Before optimizing, we check the adjoint gradient against a central finite
# difference on one degree of freedom from each class. `well_control_optimization_problem`
# returns a standalone problem object: `prob(x)` returns `(objective, gradient)`
# at a scaled control vector `x`, and `prob.x0` is the scaled base schedule.
prob = JutulDarcy.well_control_optimization_problem(copt)
obj0, grad0 = prob(prob.x0)

## Index of the first DOF of each (well, quantity) block in the control vector
dof_start = Dict{Tuple{Symbol,Symbol},Int}()
let off = 0
    for d in copt.dofs
        dof_start[(d.well, d.quantity)] = off + 1
        off += d.constant ? 1 : length(d.periods)
    end
end

function central_fd(prob, x, i; h = 1e-3)
    xp = copy(x); xp[i] += h
    xm = copy(x); xm[i] -= h
    (first(prob(xp; gradient = false)) - first(prob(xm; gradient = false))) / (2h)
end

println("gradient check at the base schedule (adjoint vs central FD):")
for (w, q) in ((:INJECT1, :rate), (:INJECT3, :rate), (:PROD1, :bhp), (:PROD3, :bhp), (:INJECT1, :bhp))
    i = dof_start[(w, q)] + 3          # period 4
    ga = grad0[i]
    gfd = central_fd(prob, prob.x0, i)
    @printf("  %-8s %-4s  adjoint % .4e   FD % .4e\n", w, q, ga, gfd)
end
# The rate and BHP *target* entries match the finite difference. The injector
# `bhp` *limit* entry is an exact zero: at the base schedule the injectors run
# at ~80 m³/day, far below their 420 bar cap, so that limit cannot influence the
# objective. It stays in the control vector and would contribute on any
# evaluation where it binds.

# ## Run the optimization
# `optimize_well_controls` drives the box-constrained quasi-Newton optimizer
# (`Jutul.unit_box_bfgs`, the same routine as MRST's `unitBoxBFGS`) with
# adjoint gradients. On the coarse model this converges in about a dozen
# iterations to roughly a 15 % NPV gain.
prm_opt = optimize_well_controls(copt; maximize = true, max_it = 15)
opt_case = copt(prm_opt)

# ## NPV over the iterations
npv_hist = copt.dopt.history.objectives
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Optimizer evaluation", ylabel = "NPV (\$100M)",
    title = @sprintf("NPV: %.3f -> %.3f  (+%.1f%%)",
        npv_hist[1], npv_hist[end], 100*(npv_hist[end] - npv_hist[1])/npv_hist[1]))
scatterlines!(ax, npv_hist)
fig

# ## Injector rate allocation
# The optimized injection rates as a fraction of the 150 m³/day upper bound, per
# injector and control period, next to the (constant) base schedule. The
# optimizer front-loads injection — near the rate ceiling early to pressurize
# and displace oil while it is cheap to produce, then tapering off as continued
# injection would only recirculate water at a cost.
prm_base = copt.dopt.parameters
rate_mat(prm) = [prm["INJECT$i"]["rate"][p] * day / 150.0 for p in 1:10, i in 1:8]

fig = Figure(size = (1100, 400))
for (j, (ttl, prm)) in enumerate(("Base schedule" => prm_base, "Optimized" => prm_opt))
    ax = Axis(fig[1, j], xlabel = "Control period", ylabel = "Injector",
        title = ttl, yticks = (1:8, string.(1:8)))
    hm = heatmap!(ax, 1:10, 1:8, rate_mat(prm), colorrange = (0, 1), colormap = :viridis)
    j == 2 && Colorbar(fig[1, 3], hm, label = "rate / 150 m³/day")
end
fig

# ## Producer BHP targets
# Each producer's optimized BHP target per control period. The producers that
# keep flowing are drawn down hard early (BHP at its 380 bar floor, maximum oil
# while the water cut is low) and eased off later (higher BHP) as water breaks
# through. The producers that water out first are instead pushed to their
# 400 bar ceiling — minimum drawdown — or shut entirely (next section).
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "Control period", ylabel = "Producer BHP target (bar)")
for (k, w) in enumerate(producers)
    scatterlines!(ax, prm_opt["PROD$k"]["bhp"] ./ bar, label = string(w))
end
hlines!(ax, [380, 400], color = :gray, linestyle = :dash)
axislegend(ax, position = :rb)
fig

# ## Which controls were active during the simulation
# For every producer we count how many report steps ran on each operating
# control. A step tagged `:disabled` is one where the well was shut: the harder
# drawdown the optimizer asks for eventually pushes a low-rate producer below
# its minimum-rate limit and it is closed for the rest of that period. On this
# coarse model the optimized schedule shuts one producer for about half the
# horizon and two others intermittently once water breaks through, while the
# fourth stays open throughout. The adjoint follows this switch — on `:disabled`
# steps the producer's BHP target correctly gets no gradient (verified against
# finite differences in the JutulDarcy test suite).
ws_base, states_base = simulate_reservoir(egg_case, info_level = -1)
ws_opt,  states_opt  = simulate_reservoir(opt_case, info_level = -1)
for w in producers
    counts = [(c, count(==(c), ws_opt.wells[w][:control])) for c in unique(ws_opt.wells[w][:control])]
    println("  $w operating controls: ", counts)
end

# ## Field production, base vs optimized
# The optimizer trades a few percent of oil for a large cut in water handling:
# on the coarse model, roughly -3 % oil, -63 % produced water and -39 % injected
# water. Water costs money both to inject and to lift, so with a 10 %/year
# discount this raises NPV by about 15 %.
mb = reservoir_measurables(egg_case.model, ws_base)
mo = reservoir_measurables(opt_case.model, ws_opt)
t_days = ws_base.time ./ day

fig = Figure(size = (1100, 400))
for (j, (key, ylab)) in enumerate((:fopt => "Field oil produced (stb)",
                                   :fwpt => "Field water produced (stb)",
                                   :fwit => "Field water injected (stb)"))
    ax = Axis(fig[1, j], xlabel = "Time (days)", ylabel = ylab)
    lines!(ax, t_days, mb[key].values ./ stb, label = "base")
    lines!(ax, t_days, mo[key].values ./ stb, label = "optimized")
    j == 1 && axislegend(ax, position = :lt)
end
fig

# ## Water saturation change
# Where the optimized schedule swept differently.
reservoir = reservoir_domain(egg_case.model)
g = physical_representation(reservoir)
sw_base = states_base[end][:Saturations][1, :]
sw_opt  = states_opt[end][:Saturations][1, :]

fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1], zreversed = true, title = "Sw(optimized) - Sw(base), final step")
plt = plot_cell_data!(ax, g, sw_opt .- sw_base, colormap = :balance, colorrange = (-0.4, 0.4))
for (w, wd) in get_model_wells(egg_case.model)
    plot_well!(ax, g, wd, top_factor = 0.8, fontsize = 8)
end
Colorbar(fig[2, 1], plt, vertical = false)
fig

# ## Differences from the MRST example
# * Single realization instead of the four-member ensemble; there is no
#   expected-value objective.
# * MRST also puts a per-control-step economic shut-in logic on each well (close
#   when the water cut exceeds 0.95, on sign change, or below 0.1 m³/day). Here
#   the only shut-in is JutulDarcy's built-in minimum-rate limit; a water-cut
#   trigger would need a custom control-logic hook.
# * MRST uses a generalized-reduced-gradient method; here we use the ported
#   `unitBoxBFGS`. Both consume adjoint gradients for the targets and for any
#   operating limit that is active during at least one simulation.
