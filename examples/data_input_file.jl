# # Simulating Eclipse/DATA input files
# The DATA format is commonly used in reservoir simulation. JutulDarcy can set
# up cases on this format and includes a fully featured grid builder for
# corner-point grids. Once a case has been set up, it uses the same types as a
# regular JutulDarcy simulation, allowing modification and use of the case in
# differentiable workflows.
#
# We begin by loading the SPE9 dataset via the GeoEnergyIO package.
using JutulDarcy, GeoEnergyIO
pth = GeoEnergyIO.test_input_file_path("SPE9", "SPE9.DATA")
# ## Set up and run a simulation
# If we do not need the case, we could also have done:
# ws, states = simulate_data_file(pth)
case = setup_case_from_data_file(pth)
ws, states = simulate_reservoir(case)
# ## Show the input data
# The input data takes the form of a Dict:
case.input_data
# We can also examine the for example RUNSPEC section, which is also represented
# as a Dict.
case.input_data["RUNSPEC"]
# ## Plot the simulation model
# These plot are interactive when run outside of the documentations.
using GLMakie
plot_reservoir(case.model, states)
# ## Plot the well responses
# We can plot the well responses (rates and pressures) in an interactive viewer.
plot_well_results(ws)
# ## Plot the field responses
# Similar to the wells, we can also plot field-wide measurables. We plot the
# field gas production rate and the average pressure as the initial selection.
plot_reservoir_measurables(case, ws, states, left = :fgpr, right = :pres)
