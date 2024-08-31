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
# ## Set up simulation model
# We set up a domain and a single injector. We pass the special :co2brine
# argument in place of the system to the reservoir model setup routine. This
# will automatically set up a compositional two-component CO2-H2O model with the
# appropriate functions for density, viscosity and miscibility.
#
# Note that this model can be run with a thermal mode by setting 
domain = reservoir_domain(mesh, permeability = 0.3Darcy, porosity = 0.3, temperature = convert_to_si(30.0, :Celsius))
Injector = setup_well(domain, (65, 1, 1), name = :Injector)
model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);
# ## Define approximate hydrostatic pressure
# The initial pressure of the water-filled domain is assumed to be at
# hydrostatic equilibrium.
nc = number_of_cells(mesh)
p0 = zeros(nc)
depth = domain[:cell_centroids][3, :]
g = Jutul.gravity_constant
@. p0 = 200bar + depth*g*1000.0

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
# ## Set up initial state
state0 = setup_reservoir_state(model,
    Pressure = p0,
    OverallMoleFractions = [1.0, 0.0],
)
# ## Simulate the schedule
# We set a maximum internal time-step of 30 days to ensure smooth convergence
# and reduce numerical diffusion.
wd, states, t = simulate_reservoir(state0, model, dt,
    parameters = parameters,
    forces = forces,
    max_timestep = 90day
)
# ## Plot the density of brine
# The density of brine depends on the CO2 concentration and gives a good
# visualization of where the mass of CO2 exists.
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
# ## Plot result in interactive viewer
plot_reservoir(model, states)