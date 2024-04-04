using Jutul, JutulDarcy, HYPRE
using GLMakie

nx = 100
nz = 50

Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)

cart_dims = (nx, 1, nz)
physical_dims = (1000.0, 1.0, 50.0)
mesh = UnstructuredMesh(CartesianMesh(cart_dims, physical_dims))

points = mesh.node_points
for (i, pt) in enumerate(points)
    x, y, z = pt
    x_u = 2*π*x/1000.0
    w = 0.2
    dz = 0.05*x + w*(30*sin(2.0*x_u) + 20*sin(5.0*x_u) + 10*sin(10.0*x_u) + 5*sin(25.0*x_u))
    points[i] = pt + [0, 0, dz]
end

boundary = Int[]
for cell in 1:number_of_cells(mesh)
    I, J, K = cell_ijk(mesh, cell)
    if I == 1 || I == nx
        push!(boundary, cell)
    end
end
##
domain = reservoir_domain(mesh, permeability = 1.0Darcy, porosity = 0.3)
Injector = setup_well(domain, (65, 1, 1), name = :Injector)
model, parameters = setup_reservoir_model(domain, :co2brine, wells = Injector);
parameters[:Reservoir][:FluidVolume][boundary] *= 1000;
##
plot_reservoir(model)
##

nstep = 25
nstep_shut = 25
dt_inject = fill(365.0day, nstep)
pv = pore_volume(model, parameters)
inj_rate = 0.01*sum(pv)/sum(dt)
rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = 700.0)
controls = Dict(:Injector => I_ctrl)
forces_inject = setup_reservoir_forces(model, control = controls)

forces_shut = setup_reservoir_forces(model)
dt_shut = fill(365.0day, nstep_shut)

dt = vcat(dt_inject, dt_shut)
forces = vcat(fill(forces_inject, nstep), fill(forces_shut, nstep_shut))

state0 = setup_reservoir_state(model,
    Pressure = 200bar,
    OverallMoleFractions = [1.0, 0.0]
)

wd, states, t = simulate_reservoir(state0, model, dt,
    parameters = parameters,
    forces = forces,
    max_timestep = 30day
)
##
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
    plot_co2!(fig, i, states[step][:PhaseMassDensities][1, :], "Brine density report step $step/$(nstep+nstep_shut)")
end

fig


##
plot_reservoir(model, states)