# # Intro to compositional flow
# This is a simple conceptual example demonstrating how to solve compositional
# flow. This example uses a two-component water-CO2 system. Note that the
# default Peng-Robinson is not accurate for this system without adjustments to
# the parameters. However, the example demonstrates the conceptual workflow for
# getting started with compositional simulation.
# ## Set up mixture
# We load the external flash package and define a two-component H2O-CO2 system.
# The constructor for each species takes in molecular weight, critical pressure,
# critical temperature, critical volume, acentric factor given as strict SI.
# This means, for instance, that molar masses are given in kg/mole and not
# g/mole or kg/kmol.
using MultiComponentFlash
h2o = MolecularProperty(0.018015268, 22.064e6, 647.096, 5.595e-05, 0.3442920843)
co2 = MolecularProperty(0.0440098, 7.3773e6, 304.1282, 9.412e-05, 0.22394)

bic = zeros(2, 2)

mixture = MultiComponentMixture([h2o, co2], A_ij = bic, names = ["H2O", "CO2"])
eos = GenericCubicEOS(mixture, PengRobinson())
# ## Set up domain and wells
using Jutul, JutulDarcy, GLMakie
nx = 50
ny = 1
nz = 20
dims = (nx, ny, nz)
g = CartesianMesh(dims, (100.0, 10.0, 10.0))
nc = number_of_cells(g)
Darcy, bar, kg, meter, Kelvin, day, sec = si_units(:darcy, :bar, :kilogram, :meter, :Kelvin, :day, :second)
K = repeat([0.1, 0.1, 0.001]*Darcy, 1, nc)
res = reservoir_domain(g, porosity = 0.3, permeability = K)
# Set up a vertical well in the first corner, perforated in top layer
prod = setup_well(g, K, [(nx, ny, 1)], name = :Producer)
# Set up an injector in the opposite corner, perforated in bottom layer
inj = setup_well(g, K, [(1, 1, nz)], name = :Injector)
# ## Define system and realize on grid
rhoLS = 844.23*kg/meter^3
rhoVS = 126.97*kg/meter^3
rhoS = [rhoLS, rhoVS]
L, V = LiquidPhase(), VaporPhase()
sys = MultiPhaseCompositionalSystemLV(eos, (L, V))
model, parameters = setup_reservoir_model(res, sys, wells = [inj, prod]);
push!(model[:Reservoir].output_variables, :Saturations)
kr = BrooksCoreyRelativePermeabilities(sys, 2.0, 0.0, 1.0)
model = replace_variables!(model, RelativePermeabilities = kr)
T0 = fill(303.15*Kelvin, nc)
parameters[:Reservoir][:Temperature] = T0
state0 = setup_reservoir_state(model, Pressure = 50*bar, OverallMoleFractions = [1.0, 0.0]);

# ## Define schedule
# 5 year (5*365.24 days) simulation period
dt0 = fill(1*day, 26)
dt1 = fill(10.0*day, 180)
dt = cat(dt0, dt1, dims = 1)
rate_target = TotalRateTarget(9.5066e-06*meter^3/sec)
I_ctrl = InjectorControl(rate_target, [0, 1], density = rhoVS)
bhp_target = BottomHolePressureTarget(50*bar)
P_ctrl = ProducerControl(bhp_target)

controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl
forces = setup_reservoir_forces(model, control = controls)
ws, states = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces);
# ## Once the simulation is done, we can plot the states
# Note that this example is intended for static publication in the
# documentation. For interactive visualization you can use functions like
# `plot_interactive` to interactively visualize the states.
z = states[end][:OverallMoleFractions][2, :]
function plot_vertical(x, t)
    data = reshape(x, (nx, nz))
    data = data[:, end:-1:1]
    fig, ax, plot = heatmap(data)
    ax.title = t
    Colorbar(fig[1, 2], plot)
    fig
end;
# ### Plot final CO2 mole fraction
plot_vertical(z, "CO2")

# ### Plot final vapor saturation 
sg = states[end][:Saturations][2, :]
plot_vertical(sg, "Vapor saturation")

# ### Plot final pressure
p = states[end][:Pressure]
plot_vertical(p./bar, "Pressure [bar]")

# ### Plot in interactive viewer
plot_reservoir(model, states, step = length(dt), key = :Saturations)
