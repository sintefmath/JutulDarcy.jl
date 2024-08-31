# # Carbon dioxoide injection in aquifer
# This example demonstrates a custom K-value compositional model for the
# injection of CO2 into a saline aquifer. The physical model for flow of CO2 is
# a realization of the description in [11th SPE Comparative Solutions
# Project](https://spe.org/en/csp/). Simulation of CO2 can be challenging, and
# we load the HYPRE package to improve performance.
using Jutul, JutulDarcy
using HYPRE
using GLMakie
nx = 100
nz = 50
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)
# ## Set up a 2D aquifer model
# We set up a Cartesian mesh that is then transformed into an unstructured mesh.
# We can then modify the coordinates to create a domain with a undulating top
# surface. CO2 will flow along the top surface and the topography of the top
# surface has a large impact on where the CO2 migrates.
cart_dims = (nx, 1, nz)
physical_dims = (1000.0, 1.0, 50.0)
mesh = UnstructuredMesh(CartesianMesh(cart_dims, physical_dims))

points = mesh.node_points
for (i, pt) in enumerate(points)
    x, y, z = pt
    x_u = 2*π*x/1000.0
    w = 0.2
    dz = 0.05*x + 0.05*abs(x - 500.0)+ w*(30*sin(2.0*x_u) + 20*sin(5.0*x_u))# + 10*sin(10.0*x_u) + 5*sin(25.0*x_u))
    # dz = 0.05*x + w*(30*sin(2.0*x_u) + 20*sin(5.0*x_u) + 10*sin(10.0*x_u) + 5*sin(25.0*x_u))
    points[i] = pt + [0, 0, dz]
end
# ## Find and plot cells intersected by a deviated injector well
# We place a single injector well. This well was unfortunately not drilled
# completely straight, so we cannot directly use `add_vertical_well` based on
# logical indices. We instead define a matrix with three columns x, y, z that
# lie on the well trajectory and use utilities from `Jutul` to find the cells
# intersected by the trajectory.
import Jutul: find_enclosing_cells, plot_mesh_edges
trajectory = [
    745.0 0.5 45;    # First point
    760.0 0.5 70;    # Second point
    810.0 0.5 100.0  # Third point
]

wc = find_enclosing_cells(mesh, trajectory)

fig, ax, plt = plot_mesh_edges(mesh, z_is_depth = true)
plot_mesh!(ax, mesh, cells = wc, transparency = true, alpha = 0.3)
lines!(ax, trajectory', color = :red)
fig

# ## Set up simulation model
# We set up a domain and a single injector. We pass the special :co2brine
# argument in place of the system to the reservoir model setup routine. This
# will automatically set up a compositional two-component CO2-H2O model with the
# appropriate functions for density, viscosity and miscibility.
#
# Note that this model by default is isothermal, but we still need to specify a
# temperature when setting up the model. This is because the properties of CO2
# strongly depend on temperature, even when thermal transport is not solved.
domain = reservoir_domain(mesh, permeability = 0.3Darcy, porosity = 0.3, temperature = convert_to_si(30.0, :Celsius))
Injector = setup_well(domain, wc, name = :Injector, simple_well = true)

model = setup_reservoir_model(domain, :co2brine, wells = Injector, extra_out = false);
# ## Customize model by adding relative permeability with hysteresis
# We define three relative permeability functions: kro(so) for the brine/liquid
# phase and krg(g) for both drainage and imbibition. Here we limit the
# hysteresis to only the non-wetting gas phase, but either combination of
# wetting or non-wetting hysteresis is supported.
#
# Note that we import a few utilities from JutulDarcy that are not exported by
# default since hysteresis falls under advanced functionality.
import JutulDarcy: table_to_relperm, add_relperm_parameters!, brooks_corey_relperm
so = range(0, 1, 10)
krog_t = so.^2
krog = PhaseRelativePermeability(so, krog_t, label = :og)

# Higher resolution for second table
sg = range(0, 1, 50)

# Evaluate Brooks-Corey to generate tables
tab_krg_drain = brooks_corey_relperm.(sg, n = 2, residual = 0.1)
tab_krg_imb = brooks_corey_relperm.(sg, n = 3, residual = 0.25)

krg_drain  = PhaseRelativePermeability(sg, tab_krg_drain, label = :g)
krg_imb  = PhaseRelativePermeability(sg, tab_krg_imb, label = :g)

fig, ax, plt = lines(sg, tab_krg_drain, label = "krg drainage")
lines!(ax, sg, tab_krg_imb, label = "krg imbibition")
lines!(ax, 1 .- so, krog_t, label = "kro")
axislegend()
fig
## Define a relative permeability variable
# JutulDarcy uses type instances to define how different variables inside the
# simulation are evaluated. The `ReservoirRelativePermeabilities` type has
# support for up to three phases with w, ow, og and g relative permeabilities
# specified as a function of their respective phases. It also supports
# saturation regions.
#
# Note: If regions are used, all drainage curves come first followed by equal
# number of imbibition curves. Since we only have a single (implicit) saturation
# region, the krg input should have two entries: One for drainage, and one for
# imbibition.
#
# We also call `add_relperm_parameters` to the model. This makes sure that when
# hysteresis is enabled, we track maximum saturation for hysteresis in each
# reservoir cell.
import JutulDarcy: KilloughHysteresis, ReservoirRelativePermeabilities
krg = (krg_drain, krg_imb) 
H_g = KilloughHysteresis() # Other options: CarlsonHysteresis, JargonHysteresis
relperm = ReservoirRelativePermeabilities(g = krg, og = krog, hysteresis_g = H_g)
replace_variables!(model, RelativePermeabilities = relperm)
add_relperm_parameters!(model)
# ## Define approximate hydrostatic pressure
# The initial pressure of the water-filled domain is assumed to be at
# hydrostatic equilibrium.
nc = number_of_cells(mesh)
p0 = zeros(nc)
depth = domain[:cell_centroids][3, :]
g = Jutul.gravity_constant
@. p0 = 200bar + depth*g*1000.0
# ## Set up initial state and parameters
state0 = setup_reservoir_state(model,
    Pressure = p0,
    OverallMoleFractions = [1.0, 0.0],
)
parameters = setup_parameters(model)

# ## Find the boundary and apply a constant pressureboundary condition
# We find cells on the left and right boundary of the model and set a constant
# pressure boundary condition to represent a bounding aquifer that retains the
# initial pressure far away from injection.

boundary = Int[]
for cell in 1:nc
    I, J, K = cell_ijk(mesh, cell)
    if I == 1 || I == nx
        push!(boundary, cell)
    end
end
bc = flow_boundary_condition(boundary, domain, p0[boundary], fractional_flow = [1.0, 0.0])

# ## Plot the model
plot_reservoir(model)
# ## Set up schedule
# We set up 25 years of injection and 25 years of migration where the well is
# shut. The density of the injector is set to 900 kg/m^3, which is roughly
# the density of CO2 at in-situ conditions.
nstep = 25
nstep_shut = 25
dt_inject = fill(365.0day, nstep)
pv = pore_volume(model, parameters)
inj_rate = 0.05*sum(pv)/sum(dt_inject)

rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0],
    density = 900.0,
)
# Set up forces for use in injection
controls = Dict(:Injector => I_ctrl)
forces_inject = setup_reservoir_forces(model, control = controls, bc = bc)
# Forces with shut wells
forces_shut = setup_reservoir_forces(model, bc = bc)
dt_shut = fill(365.0day, nstep_shut);
# Combine the report steps and forces into vectors of equal length
dt = vcat(dt_inject, dt_shut)
forces = vcat(
    fill(forces_inject, nstep),
    fill(forces_shut, nstep_shut)
);
# ## Add some more outputs for plotting
rmodel = reservoir_model(model)
push!(rmodel.output_variables, :RelativePermeabilities)
push!(rmodel.output_variables, :PhaseViscosities)
# ## Simulate the schedule
# We set a maximum internal time-step of 30 days to ensure smooth convergence
# and reduce numerical diffusion.
wd, states, t = simulate_reservoir(state0, model, dt,
    parameters = parameters,
    forces = forces,
    max_timestep = 90day
)
# ## Plot the CO2 mole fraction
# We plot log10 of the CO2 mole fraction. We use log10 to account for the fact
# that the mole fraction in cells made up of only the aqueous phase is much
# smaller than that of cells with only the gaseous phase, where there is almost
# just CO2.
using GLMakie
function plot_co2!(fig, ix, x, title = "")
    ax = Axis3(fig[ix, 1],
        zreversed = true,
        azimuth = -0.51π,
        elevation = 0.05,
        aspect = (1.0, 1.0, 0.3),
        title = title)
    plt = plot_cell_data!(ax, mesh, x, colormap = :seaborn_icefire_gradient)
    Colorbar(fig[ix, 2], plt)
end
fig = Figure(size = (900, 1200))
for (i, step) in enumerate([1, 5, nstep, nstep+nstep_shut])
    plot_co2!(fig, i, log10.(states[step][:OverallMoleFractions][2, :]), "log10 of CO2 mole fraction at report step $step/$(nstep+nstep_shut)")
end
fig
# ## Plot all relative permeabilities for all time-steps
# We can plot all relative permeability evaluations. This both verifies that the
# hysteresis model is active, but also gives an indication to how many cells are
# exhibiting imbibition during the simulation.
kro_val = Float64[]
krg_val = Float64[]
sg_val = Float64[]
for state in states
    kr_state = state[:RelativePermeabilities]
    s_state = state[:Saturations]
    for c in 1:nc
        push!(kro_val, kr_state[1, c])
        push!(krg_val, kr_state[2, c])
        push!(sg_val, s_state[2, c])
    end
end

fig = Figure()
ax = Axis(fig[1, 1], title = "Relative permeability during simulation")
fig, ax, plt = scatter(sg_val, kro_val, label = "kro", alpha = 0.3)
scatter!(ax, sg_val, krg_val, label = "krg", alpha = 0.3)
axislegend()
fig
# ## Plot result in interactive viewer
# If you have interactive plotting available, you can explore the results
# yourself.
plot_reservoir(model, states)
