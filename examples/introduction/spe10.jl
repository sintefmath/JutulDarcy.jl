# # SPE10, model 2
# <tags: Introduction, InputFile>
# The SPE10 benchmark case is a standard benchmark case for reservoir
# simulators. The model is described in detail in the [SPE10 benchmark
# page](https://www.spe.org/web/csp/datasets/set02.htm)
#
# This example demonstrates routines to set up the SPE10, model 2, case using
# premade functions in the `JutulDarcy.SPE10` module. This model is often just
# called SPE10, since the first model of the benchmark is very small.
# ## Set up the reservoir
# We can set up the reservoir itself to have a look at the static properties.
using JutulDarcy, GLMakie, HYPRE
reservoir = JutulDarcy.SPE10.setup_reservoir()
plot_reservoir(reservoir, key = :porosity)
# ## Set up and run a simulation for the first layer
# We can set up and run a full simulation for the first layer of the model. As
# we want the example to run quickly, we just pick the top layer. The function
# is set up to scale the default well rates to the smaller model, so you can
# adjust the number of layers without changing anything else.
case = JutulDarcy.SPE10.setup_case(layers = 1:1)
# ### Plot the porosity
# Note that some cells are removed due to very low porosity. The setup function
# has additional options to control this behavior if you would rather limit the minimum
# porosity than removing cells.
plot_reservoir(case.model, key = :porosity)
# ### Run the simulation
ws, states = simulate_reservoir(case)
# ### Plot the final saturation
plot_reservoir(case.model, states, key = :Saturations, step = length(states))
# ## Show the last layer
# The first 35 layers correspond to the Tarbert formation and the latter 50
# layers are the upper ness formation. We plot one of the last layers, showing a
# channelized, fluvial structure.
reservoir_ness = JutulDarcy.SPE10.setup_reservoir(layers = 60)
plot_reservoir(reservoir_ness, key = :porosity)
# ## Coarsen the model
# The full model is quite large, with 1.1 million cells. We can coarsen the
# model to get a smaller model that is faster to simulate. The original model
# was intended as a benchmark for upscaling methods, even though such models are
# now quite easy to simulate when using parallel computing on CPU/GPU.
#
# The default coarsening is quite simple (harmonic averages and sums), but
# quickly sets up a coarse model that can be used for testing.
case_fine = JutulDarcy.SPE10.setup_case()
case_coarse = coarsen_reservoir_case(case_fine, (10, 22, 8))
# ### Simulate the coarse model
# The now quite coarse model should run in less than a second.
ws_coarse, states_coarse = simulate_reservoir(case_coarse);
plot_reservoir(case_coarse.model, states_coarse, key = :Saturations, step = length(states_coarse))
# ## Conclusion
# The SPE10, model 2, case is a standard benchmark case for reservoir
# simulators. JutulDarcy includes premade functions to set up and run this
# model, including wells, PVT, and the schedule from the original case. The
# model can be run as-is, or coarsened to a smaller size for testing and
# development.
