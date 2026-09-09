# # Well control optimization for NPV (Egg model)
# <tags: Immiscible, InputFile, ModelReduction, Optimization, Advanced>
#
# We optimize the well controls of the (coarsened) Egg model over a sequence of
# control periods to maximize net present value. Compared with
# [the rate optimization example](rate_optimization.md) — one injector rate per
# report step with a hand-built objective — this one uses the general
# [`setup_well_control_optimization`](@ref) helper and an *extended control
# vector* that mixes well *targets* (injector rates, producer BHPs) with
# operating *limits* (an injector BHP cap) in the same optimization. The adjoint
# follows whichever control is active during the simulation, so a limit
# contributes to the gradient exactly on the steps where it binds.

using Jutul, JutulDarcy, GLMakie, GeoEnergyIO, LinearAlgebra, Statistics, Printf

# ## Load the Egg model
# First three years of the Egg schedule, coarsened to 20x20x3 for speed. Use
# `egg_case = fine_case[1:n5]` to optimize on the full grid.
egg_dir = GeoEnergyIO.test_input_file_path("EGG")
fine_case = setup_case_from_data_file(joinpath(egg_dir, "EGG.DATA"))
n5 = count(<=(3*si_unit(:year)), cumsum(fine_case.dt))
egg_case = coarsen_reservoir_case(fine_case[1:n5], (20, 20, 3), method = :ijk)
nstep = length(egg_case.dt)

# ## The base schedule
# Eight water injectors under total-rate control (79.5 m³/day, 420 bar BHP cap)
# and four producers under BHP control (395 bar), one control for the whole
# horizon.
injectors = [Symbol("INJECT$i") for i in 1:8]
producers = [Symbol("PROD$i")   for i in 1:4]

day  = si_unit(:day)
bar  = si_unit(:bar)
stb  = si_unit(:stb)

# ## Control periods
# Split the horizon into `n_control_periods` groups of report steps of roughly
# equal length. Each well quantity gets an independent value in each period.
n_control_periods = 6
edges = round.(Int, range(0, nstep, length = n_control_periods+1))
t = cumsum(egg_case.dt)
step_to_control = ceil.(Int, n_control_periods * t ./ t[end])
edges = [0; findall(diff(step_to_control) .!= 0); nstep]
control_periods = [(edges[i] + 1):edges[i + 1] for i in 1:n_control_periods]

# ## The extended control vector
# `setup_well_control_optimization` takes a list of degree-of-freedom specs. A
# quantity that matches the well's base control target is optimized as the
# *target*; any other quantity becomes an operating *limit* written into
# `forces[:Facility].limits`. Bounds are absolute (`abs_min`/`abs_max`) and/or
# relative to the base value (`rel_min`/`rel_max`).
#
# | quantity | role | bounds |
# |---|---|---|
# | injector `rate` | target | `[0, 150]` m³/day |
# | injector `bhp` | limit (upper) | `[400, 420]` bar |
# | producer `bhp` | target | `[380, 400]` bar |
#
# Injectors are rate-controlled and producers BHP-controlled, so the injector
# `bhp` entry is an operating limit: it contributes to the gradient only on the
# steps of the evaluations where the well actually sits on the cap, and is inert
# otherwise.
controls = NamedTuple[]
for w in injectors
    push!(controls, (well = w, quantity = :rate, abs_min = 0.0/day,  abs_max = 150.0/day))
    push!(controls, (well = w, quantity = :bhp,  abs_min = 400*bar,  abs_max = 420*bar))
end
for w in producers
    push!(controls, (well = w, quantity = :bhp,  abs_min = 380*bar,  abs_max = 400*bar))
end

# ## The NPV objective
# [`JutulDarcy.npv_objective`](@ref): oil worth 50 \$/stb, injected and produced
# water each cost 2 \$/stb, 10 %/year discount rate. `scale` brings the objective
# to order unity for the optimizer.
npv_arg = (
    oil_price     = 50.0,
    water_price   = -2.0,   # produced water: a cost
    water_cost    =  2.0,   # injected water: a cost
    oil_cost      =  0.0,
    discount_rate =  0.1,
    scale         =  1e8,
)

copt = setup_well_control_optimization(egg_case, control_periods, controls;
    objective = :npv, npv_arg = npv_arg)

# ## Run the optimization
# `optimize_well_controls` drives a box-constrained quasi-Newton optimizer with
# adjoint gradients (sparse differentiation over the control vector by default).
prm_opt = optimize_well_controls(copt; maximize = true, max_it = 15,
    simulator_arg = (rtol = 1e-5, tol_cnv = 1e-5))
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
# Optimized injection rate as a fraction of the 150 m³/day upper bound, per
# injector and control period, next to the constant base schedule.
prm_base = copt.dopt.parameters
rate_mat(prm) = [prm["INJECT$i"]["rate"][p] * day / 150.0 for p in 1:n_control_periods, i in 1:8]

fig = Figure(size = (1100, 400))
for (j, (ttl, prm)) in enumerate(("Base schedule" => prm_base, "Optimized" => prm_opt))
    local ax = Axis(fig[1, j], xlabel = "Control period", ylabel = "Injector",
        title = ttl, yticks = (1:8, string.(1:8)))
    hm = heatmap!(ax, 1:n_control_periods, 1:8, rate_mat(prm), colorrange = (0, 1), colormap = :viridis)
    j == 2 && Colorbar(fig[1, 3], hm, label = "rate / 150 m³/day")
end
fig

# ## Producer BHP targets
# Each producer's optimized BHP target per control period, between the 380 and
# 400 bar bounds — typically a hard drawdown early and easing off as water breaks
# through.
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "Control period", ylabel = "Producer BHP target (bar)")
for (k, w) in enumerate(producers)
    scatterlines!(ax, prm_opt["PROD$k"]["bhp"] ./ bar, label = string(w))
end
hlines!(ax, [380, 400], color = :gray, linestyle = :dash)
axislegend(ax, position = :rb)
fig

# ## Which controls were active during the simulation
# Report-step count per operating control for each well. `:bhp` on an injector
# means it hit its BHP cap and switched off rate control; `:disabled` means the
# well was shut after its rate dropped below the minimum-rate limit. The adjoint
# follows these switches — see `test/well_control_gradients.jl`.
ws_base, states_base = simulate_reservoir(egg_case, info_level = -1)
ws_opt,  states_opt  = simulate_reservoir(opt_case, info_level = -1)
for w in [injectors; producers]
    counts = [(c, count(==(c), ws_opt.wells[w][:control])) for c in unique(ws_opt.wells[w][:control])]
    println("  $w operating controls: ", counts)
end

# ## Injectors: simulation vs optimized schedule vs bounds
# For a few injectors, over time: the simulated rate/BHP, the optimized schedule
# value (a staircase — the rate *target* and the BHP *limit*), and the
# `[abs_min, abs_max]` box. Where the simulated rate falls below the target the
# well has hit its BHP cap; the simulated BHP then sits on that cap.
sel_injectors = [:INJECT1, :INJECT6, :INJECT8]

tt = ws_opt.time ./ day
function period_span(p)                       # (start, end) time of control period p, in days
    s = control_periods[p]
    (first(s) == 1 ? 0.0 : tt[first(s) - 1], tt[last(s)])
end
dof_box(well, quantity) =                      # [abs_min, abs_max] for a (well, quantity) DOF
    let d = only(filter(d -> d.well == well && d.quantity == quantity, copt.dofs))
        (d.abs_min, d.abs_max)
    end

# ### Water (surface) injection rate and BHP
fig = Figure(size = (360 * length(sel_injectors), 380))
for (j, w) in enumerate(sel_injectors)
    local ax = Axis(fig[1, j], xlabel = "Time (days)", ylabel = "Injection rate (m³/day)",
        title = string(w))
    lo, hi = dof_box(w, :rate)
    hlines!(ax, [lo, hi] .* day, color = :gray, linestyle = :dash, label = "optimization bounds")
    for p in 1:n_control_periods
        t0, t1 = period_span(p)
        v = prm_opt[string(w)]["rate"][p] * day
        lines!(ax, [t0, t1], [v, v], color = :red, linewidth = 3, linestyle = :dash,
            label = p == 1 ? "rate target (optimized)" : nothing)
    end
    lines!(ax, tt, abs.(ws_opt.wells[w][:rate]) .* day, color = Cycled(1), linewidth = 2,
        label = "simulated")
    j == 1 && axislegend(ax, position = :lb, framevisible = false)

    ax = Axis(fig[2, j], xlabel = "Time (days)", ylabel = "BHP (bar)", title = string(w))
    lo, hi = dof_box(w, :bhp)
    hlines!(ax, [lo, hi] ./ bar, color = :gray, linestyle = :dash, label = "optimization bounds")
    for p in 1:n_control_periods
        t0, t1 = period_span(p)
        v = prm_opt[string(w)]["bhp"][p] / bar
        lines!(ax, [t0, t1], [v, v], color = :red, linewidth = 3, linestyle = :dash,
            label = p == 1 ? "BHP limit (optimized)" : nothing)
    end
    lines!(ax, tt, ws_opt.wells[w][:bhp] ./ bar, color = Cycled(1), linewidth = 2,
        label = "simulated")
    j == 1 && axislegend(ax, position = :lb, framevisible = false)
end
fig

# ## Field production, base vs optimized
# The optimizer trades a little oil for a large cut in produced and injected
# water — water costs money both to lift and to inject.
mb = reservoir_measurables(egg_case.model, ws_base)
mo = reservoir_measurables(opt_case.model, ws_opt)
t_days = ws_base.time ./ day

fig = Figure(size = (1100, 400))
for (j, (key, ylab)) in enumerate((:fopt => "Field oil produced (stb)",
                                   :fwpt => "Field water produced (stb)",
                                   :fwit => "Field water injected (stb)"))
    local ax = Axis(fig[1, j], xlabel = "Time (days)", ylabel = ylab)
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

# ## Notes
# * The only shut-in mechanism here is JutulDarcy's built-in minimum-rate limit;
#   an economic trigger (e.g. shut on water cut) would need a custom
#   control-logic hook.
# * A single realization: no ensemble / expected-value objective.
