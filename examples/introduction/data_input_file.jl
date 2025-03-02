# # Simulating Eclipse/DATA input files
# The DATA format is commonly used in reservoir simulation. JutulDarcy can set
# up cases on this format and includes a fully featured grid builder for
# corner-point grids. Once a case has been set up, it uses the same types as a
# regular JutulDarcy simulation, allowing modification and use of the case in
# differentiable workflows.
#
# We begin by loading the SPE9 dataset via the GeoEnergyIO package. This package
# includes a set of open datasets that can be used for testing and benchmarking.
# The SPE9 dataset is a 3D model with a corner-point grid and a set of wells
# produced by the Society of Petroleum Engineers. The specific version of the
# file included here is taken from the [OPM
# tests](https://github.com/OPM/opm-tests) repository.
using JutulDarcy, GeoEnergyIO
pth = GeoEnergyIO.test_input_file_path("SPE9", "SPE9.DATA");
# ## Set up and run a simulation
# We have supressed the output of the simulation to avoid cluttering the
# documentation, but we can set the `info_level` to a higher value to see the
# output.
#
# If we do not need the case, we could also have simulated by passing the path:
# `ws, states = simulate_data_file(pth)`
case = setup_case_from_data_file(pth)
ws, states = simulate_reservoir(case);
# ## Show the input data
# The input data takes the form of a Dict:
case.input_data
# We can also examine the for example RUNSPEC section, which is also represented
# as a Dict.
case.input_data["RUNSPEC"]
# ## Plot the simulation model
# These plot are normally interactive, but if you are reading the published
# online documentation static screenshots will be inserted instead.
using GLMakie
plot_reservoir(case.model, states)
# ## Plot the well responses
# We can plot the well responses (rates and pressures) in an interactive viewer.
# Multiple wells can be plotted simultaneously, with options to select which
# units are to be used for plotting.
plot_well_results(ws)
# ## Plot the field responses
# Similar to the wells, we can also plot field-wide measurables. We plot the
# field gas production rate and the average pressure as the initial selection.
# If you are running this case interactively you can select which measurables to
# plot.
#
# We observe that the field pressure steadily decreases over time, as a result
# of the gas production. The drop in pressure is not uniform, as during the
# period where little gas is produced, the decrease in field pressure is slower.
plot_reservoir_measurables(case, ws, states, left = :fgpr, right = :pres)
