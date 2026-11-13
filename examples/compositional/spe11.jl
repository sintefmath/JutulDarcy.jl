# # SPE11 benchmark models
# <tags: Introduction, Compositional>
# The [11th SPE Comparative Solutions Project](https://www.spe.org/en/csp/) is a
# benchmark on simulating CO2 storage. See the published papers for more information:
# - [Benchmark setup](https://onepetro.org/SJ/article/doi/10.2118/218015-PA/540636)
# - [Simulation results](https://doi.org/10.1016/j.ijggc.2025.104519)
#
# SPE11 is a set of three benchmark models. Unlike previous benchmarks, the grid
# resolution is not set, but rather a part of the modelling exercise. JutulDarcy
# provides scripts to set up these benchmark cases, with options to customize
# grid resolution and adjust the parameters.
#
# This script demonstrates the three cases, at very coarse resolutions.
using JutulDarcy, Jutul, GLMakie
# ## Case A
# Case A is laboratory-scale. This model is generally quite difficult to
# simulate with higher resolutions due to the near-atmospheric conditions. We
# simulate a very coarse version here.
nx = 20
ny = 20
case_a = JutulDarcy.CSP11.setup_case_a((nx, 1, ny))

plot_reservoir(case_a, key = :satnum, fancy = false)
# ### Simulate Case A
# We use ILU(0) and a slightly looser tolerance for this model, as CPR can
# produce failures near 1 atmosphere pressure.
ws_a, states_a = simulate_reservoir(case, tol_cnv = 0.01, precond = :ilu0)
plot_reservoir(case_a, states_a, key = :Saturations, step = length(case_a.dt), fancy = false)
# ## Case B
# Case B is a 2D version at reservoir conditions that also incorporates thermal effects.
case_b = JutulDarcy.CSP11.setup_case_b((50, 1, 50), thermal = true, dz_max = 0.1)
ws_b, states_b = simulate_reservoir(case_b, tol_cnv = 0.01)
plot_reservoir(case_b, states_b, key = :Saturations, step = length(case_b.dt), fancy = false)
# ## Case C
# Case C is a 3D version at reservoir conditions.
case = JutulDarcy.CSP11.setup_case_c((30, 30, 30), thermal = true, dz_max = 0.1, nstep_initialization = 0)
ws_c, states_c = simulate_reservoir(case, tol_cnv = 0.01)
##
msh = reservoir_mesh(case)
z_co2 = states_c[end][:OverallMoleFractions][2, :]
fig, ax, plt = plot_cell_data(msh, z_co2, cells = z_co2 .> 0.01)
for (k, w) in get_model_wells(case)
    plot_well!(ax, msh, w)
end
plot_mesh_edges!(ax, msh)
fig
