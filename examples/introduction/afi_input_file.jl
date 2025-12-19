# # Set up a .AFI file for simulation
# <tags: Introduction, InputFile>
#
# JutulDarcy has partial support for setting up and reading the .AFI file
# format. This example demonstrates how to load AFI files.
# ## Load the OLYMPUS model
# We first load the Olympus model from a RESQML-based AFI file.
using Jutul, JutulDarcy, GeoEnergyIO
fn = GeoEnergyIO.test_input_file_path("OLYMPUS_25_AFI_RESQML", "OLYMPUS_25.afi")
case_olympus = setup_case_from_afi(fn);
# ## Reading files with a pre-defined reservoir
# Note that as the AFI support requires either inline mesh definitions or
# RESQML, and as such, not all files can be read directly. However, if we have
# already set up a reservoir model (for example from a DATA file), we can reuse
# the mesh and mesh properties from that model when setting up the AFI case. We
# set up SPE9 with a pre-defined reservoir model, bypassing the need for GSG
# support.
spe9_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE9")
case = setup_case_from_data_file(joinpath(spe9_dir, "SPE9.DATA"))
reservoir = reservoir_domain(case)
fn = GeoEnergyIO.test_input_file_path("SPE9_AFI_GSG", "SPE9_clean_split.afi")
case_ix = setup_case_from_afi(fn, reservoir = reservoir);
