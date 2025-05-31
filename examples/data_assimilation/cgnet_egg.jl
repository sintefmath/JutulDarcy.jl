# # History matching a coarse model - CGNet
# This example demonstrates how to calibrate a coarse model against results from
# a fine model. We do this by optimizing the parameters of the coarse model to
# match the well curves. This is a implementation of the method described in
# [cgnet1](@cite). This also serves as a demonstration of how to use the
# simulator for history matching, as the fine model results can stand in for
# real field observations.
#
# ## Load and simulate Egg base case
# We take a subset of the first 60 steps (1350 days) since not much happens after
# that in terms of well behavior.
using Jutul, JutulDarcy, HYPRE, GeoEnergyIO, GLMakie
import LBFGSB as lb

egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG")
data_pth = joinpath(egg_dir, "EGG.DATA")

fine_case = setup_case_from_data_file(data_pth)
fine_case = fine_case[1:60]
simulated_fine = simulate_reservoir(fine_case)
plot_reservoir(fine_case, simulated_fine.states, key = :Saturations, step = 60)
# ## Create initial coarse model and simulate
coarse_case = JutulDarcy.coarsen_reservoir_case(fine_case, (25, 25, 5), method = :ijk)
simulated_coarse = simulate_reservoir(coarse_case)
plot_reservoir(coarse_case, simulated_coarse.states, key = :Saturations, step = 60)
# ## Setup optimization
# We set up the optimization problem by defining the objective function as a sum
# of squared mismatches for all well observations, for all time-steps. We also
# define limits for the parameters, and set up the optimization problem.
#
# We also limit the number of function evaluations since this example runs as a
# part of continuous integration and we want to keep the runtime short.
function setup_optimization_cgnet(case_c, case_f, result_f)
    states_f = result_f.result.states
    wells_results, = result_f
    model_c = case_c.model
    state0_c = setup_state(model_c, case_c.state0);
    param_c = setup_parameters(model_c)
    forces_c = case_c.forces
    dt = case_c.dt
    model_f = case_f.model

    bhp = JutulDarcy.BottomHolePressureTarget(1.0)
    wells = collect(keys(JutulDarcy.get_model_wells(case_f)))

    day = si_unit(:day)
    wrat_scale = (1/150)*day
    orat_scale = (1/80)*day
    grat_scale = (1/1000)*day

    w = []
    matches = []
    signs = []
    sys = reservoir_model(model_f).system
    wrat = SurfaceWaterRateTarget(-1.0)
    orat = SurfaceOilRateTarget(-1.0)
    grat = SurfaceGasRateTarget(-1.0)

    push!(matches, bhp)
    push!(w, 1.0/si_unit(:bar))
    push!(signs, 1)

    for phase in JutulDarcy.get_phases(sys)
        if phase == LiquidPhase()
            push!(matches, orat)
            push!(w, orat_scale)
            push!(signs, -1)

        elseif phase == VaporPhase()
            push!(matches, grat)
            push!(w, grat_scale)
            push!(signs, -1)
        else
            @assert phase == AqueousPhase()
            push!(matches, wrat)
            push!(w, wrat_scale)
            push!(signs, -1)
        end
    end
    signs = zeros(Int, length(signs))
    o_scale = 1.0/(sum(dt)*length(wells))
    G = (model_c, state_c, dt, step_info, forces) -> well_mismatch(
        matches,
        wells,
        model_f,
        states_f,
        model_c,
        state_c,
        dt,
        step_info,
        forces,
        weights = w,
        scale = o_scale,
        signs = signs
    )

    @assert Jutul.evaluate_objective(G, model_f, states_f, dt, case_f.forces) == 0.0
    ##
    cfg = optimization_config(model_c, param_c,
        use_scaling = true,
        rel_min = 0.001,
        rel_max = 1000
    )
    for (k, v) in cfg
        for (ki, vi) in v
            if ki == :FluidVolume
                vi[:active] = k == :Reservoir
            end
            if ki == :ConnateWater
                vi[:active] = false
            end
            if ki in [:TwoPointGravityDifference, :PhaseViscosities, :PerforationGravityDifference]
                vi[:active] = false
            end
            if ki in [:WellIndices, :Transmissibilities]
                vi[:active] = true
                vi[:abs_min] = 0.0
                vi[:abs_max] = 1e-6
            end

        end
    end
    opt_setup = setup_parameter_optimization(model_c, state0_c, param_c, dt, forces_c, G, cfg);
    x0 = opt_setup.x0
    F0 = opt_setup.F!(x0)
    dF0 = opt_setup.dF!(similar(x0), x0)
    println("Initial objective: $F0, gradient norm $(sum(abs, dF0))")
    return opt_setup
end

# ## Define the optimization loop
# JutulDarcy can use any optimization package that can work with gradients and
# limits, here we use the `LBFGSB` package.
function optimize_cgnet(opt_setup)
    lower = opt_setup.limits.min
    upper = opt_setup.limits.max
    x0 = opt_setup.x0
    n = length(x0)
    setup = Dict(:lower => lower, :upper => upper, :x0 => copy(x0))

    prt = 1
    f! = (x) -> opt_setup.F_and_dF!(NaN, nothing, x)
    g! = (dFdx, x) -> opt_setup.F_and_dF!(NaN, dFdx, x)
    results, final_x = lb.lbfgsb(f!, g!, x0, lb=lower, ub=upper,
        iprint = prt,
        factr = 1e12,
        maxfun = 20,
        maxiter = 20,
        m = 20
    )
    return (final_x, results, setup)
end
# ## Run the optimization
opt_setup = setup_optimization_cgnet(coarse_case, fine_case, simulated_fine);
final_x, results, setup = optimize_cgnet(opt_setup);
# ### Transfer the results back
# The optimization is generic and works on a single long vector that represents
# all our active parameters. We can devectorize this vector back into the nested
# representation used by the model itself and simulate.
tuned_case = deepcopy(opt_setup.data[:case])
model_c = coarse_case.model
model_f = fine_case.model
param_c = tuned_case.parameters
data = opt_setup.data
devectorize_variables!(param_c, model_c, final_x, data[:mapper], config = data[:config])


simulated_tuned = simulate_reservoir(tuned_case);
# ### Plot the results interactively
using GLMakie

wells_f, = simulated_fine
wells_c, = simulated_coarse
wells_t, states_t, time = simulated_tuned

plot_well_results([wells_f, wells_c, wells_t], time, names = ["Fine", "CGNet-initial", "CGNet-tuned"])
# ### Create a function to compare individual wells
# We next compare individual wells to see how the optimization has affected the
# match between the coarse scale and fine scale. As we can see, we have
# reasonably good match between the original model with about 18,000 cells and
# the coarse model with about 3000 cells. Even better match could be possible by
# adding more coarse blocks, or also optimizing for example the relative
# permeability parameters for the coarse model.
#
# We plot the water cut and total rate for the production wells, and the bottom
# hole pressure for the injection wells.
function plot_tuned_well(k, kw; lposition = :lt)
    fig = Figure()
    ax = Axis(fig[1, 1], title = "$k", xlabel = "days", ylabel = "$kw")
    t = wells_f.time./si_unit(:day)
    if kw == :wcut
        F = x -> x[k, :wrat]./x[k, :lrat]
    else
        F = x -> abs.(x[k, kw])
    end

    lines!(ax, t, F(wells_f), label = "Fine-scale")
    lines!(ax, t, F(wells_c), label = "Initial coarse")
    lines!(ax, t, F(wells_t), label = "Tuned coarse")
    axislegend(position = lposition)
    fig
end
# ### Plot PROD1 water cut
plot_tuned_well(:PROD1, :wcut)
# ### Plot PROD2 water cut
plot_tuned_well(:PROD2, :wcut)
# ### Plot PROD4 water cut
plot_tuned_well(:PROD4, :wcut)
# ### Plot PROD1 total rate
plot_tuned_well(:PROD1, :rate)
# ### Plot PROD2 total rate
plot_tuned_well(:PROD2, :rate)
# ### Plot PROD4 total rate
plot_tuned_well(:PROD4, :rate)
# ### Plot INJECT1 bhp
plot_tuned_well(:INJECT1, :bhp, lposition = :rt)
# ### Plot INJECT4 bhp
plot_tuned_well(:INJECT4, :bhp, lposition = :rt)

# ## Plot the objective evaluations during optimization
fig = Figure()
ys = log10
is = x -> x
ax1 = Axis(fig[1, 1], yscale = ys, title = "Objective evaluations", xlabel = "Iterations", ylabel = "Objective")
plot!(ax1, opt_setup[:data][:obj_hist][2:end])
fig
# ## Plot the evoluation of scaled parameters
# We show the difference between the initial and final values of the scaled
# parameters, as well as the lower bound.
#
# JutulDarcy maps the parameters to a single vector for optimization with values
# that are approximately in the box limit range (0, 1). This is convenient for
# optimizers, but can also be useful when plotting the parameters, even if the
# units are not preserved in this visualization, only the magnitude.
fig = Figure(size = (800, 600))
ax1 = Axis(fig[1, 1], title = "Scaled parameters", ylabel = "Scaled value")
scatter!(ax1, setup[:x0], label = "Initial X")
scatter!(ax1, final_x, label = "Final X", markersize = 5)
lines!(ax1, setup[:lower], label = "Lower bound")
axislegend()

trans = data[:mapper][:Reservoir][:Transmissibilities]

function plot_bracket(v, k)
    start = v.offset_x+1
    stop = v.offset_x+v.n_x
    y0 = setup[:lower][start]
    y1 = setup[:lower][stop]
    bracket!(ax1, start, y0, stop, y1,
    text = "$k", offset = 1, orientation = :down)
end

for (k, v) in pairs(data[:mapper][:Reservoir])
    plot_bracket(v, k)
end
ylims!(ax1, (-0.2*maximum(final_x), nothing))
fig

