# # Adjoint gradients for the SPE1 model
# This is a brief example of how to set up a case and compute adjoint gradients
# for the SPE1 model in a few different ways. These techniques are in general
# applicable to DATA-type models.
# ## Load the model
using Jutul, JutulDarcy, GeoEnergyIO, GLMakie
data_pth = joinpath(GeoEnergyIO.test_input_file_path("SPE1"), "SPE1.DATA")
data = parse_data_file(data_pth);
case = setup_case_from_data_file(data);
# ## Set up a function to set up the case with custom porosity
# We create a setup function that takes in a parameter dictionary `prm` and
# returns a case with the porosity set to the value in `prm["poro"]`. This is a
# convenient way to customize a case prior to setup, allowing direct interaction
# with keywords like PERMX, PORO, etc.
function F(prm, step_info = missing)
    data_c = deepcopy(data)
    data_c["GRID"]["PORO"] = fill(prm["poro"], size(data_c["GRID"]["PORO"]))
    case = setup_case_from_data_file(data_c)
    return case
end
# ## Define and simulate truth case (poro from input)
x_truth = only(unique(data["GRID"]["PORO"]))
prm_truth = Dict("poro" => x_truth)
case_truth = F(prm_truth)
ws, states = simulate_reservoir(case_truth)
# ## Define pressure difference as objective function
function pdiff(p, p0)
    v = 0.0
    for i in eachindex(p)
        v += (p[i] - p0[i])^2
    end
    return v
end

step_times = cumsum(case.dt)
function mismatch_objective(m, s, dt, step_info, forces)
    t = step_info[:time] + dt
    step = findmin(x -> abs(x - t), step_times)[2]
    p = s[:Reservoir][:Pressure]
    v = pdiff(p, states[step][:Pressure])
    return dt*(v/(si_unit(:bar)*100)^2)
end
# ## Create a perturbed initial guess and optimize
# We create a perturbed initial guess for the porosity and optimize the case
# using the mismatch objective function defined above. The optimization will
# adjust the porosity to minimize the mismatch between the simulated pressure
# and the truth case, and recover the original porosity value of 0.3.
prm = Dict("poro" => x_truth .+ 0.25)
dprm = setup_reservoir_dict_optimization(prm, F)
free_optimization_parameter!(dprm, "poro", abs_max = 1.0, abs_min = 0.1)
prm_opt = optimize_reservoir(dprm, mismatch_objective);
dprm
# ## Plot the optimization history
using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "LBFGS iteration", ylabel = "Objective function", yscale = log10)
scatter!(ax, dprm.history.val)
fig
# ## Compute sensitivities outside the optimization
# We can also compute the sensitivities outside the optimization process. As our
# previous setup funciton only has a single parameter (the porosity), we instead
# switch to the `setup_reservoir_dict_optimization` function that can set up a
# set of "typical" tunable parameters for any reservoir model. This saves us the
# hassle of writing this function ourselves when we want to optimize e.g.
# permeability, porosity and well indices.
dprm_case = setup_reservoir_dict_optimization(case)
free_optimization_parameters!(dprm_case)
dprm_grad = parameters_gradient_reservoir(dprm_case, mismatch_objective);
# ## Plot the gradient of the mismatch objective with respect to the porosity
# We see, as expected, that the gradient is largest in magnitude around the
# wells and near the front of the displacement.
m = physical_representation(reservoir_domain(case.model))
plot_cell_data(m, dprm_grad[:model][:porosity])
# ## Plot the sensitivities in the interactive viewer
# If you are running the example yourself, you can now explore the sensitivities
# in the interactive viewer. This is useful for understanding how the model
# responds to changes in the parameters.
plot_reservoir(case.model, dprm_grad[:model])
