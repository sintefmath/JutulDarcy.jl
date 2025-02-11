# # Borehole Thermal Energy Storage (BTES)
# This script demonstrates how to model a borehole thermal energy storage (BTES)
# using JutulDarcy. We will set up a BTES system with 50 wells, and simulate 10
# year of operation with 6 months of charging followed by 6 months of
# discharging.
using Jutul, JutulDarcy
using HYPRE
using GLMakie
using Gmsh

atm, kilogram, meter, Kelvin, joule, watt, litre, year, second, darcy =
    si_units(:atm, :kilogram, :meter, :Kelvin, :joule, :watt, :litre, :year, :second, :darcy)

# ## Define BTES well coordinates
# We place the BTES wells in a pattern inspired by sunflower seeds, ensuring
# that they are uniformy distributed over the domain.
num_wells = 50
radius = 15.0meter
φ = (1 + sqrt(5)) / 2 # Golden ratio
Δθ = 2 * π / φ^2
p = (k) -> radius .* sqrt(k / num_wells) .* (cos(k * Δθ), sin(k * Δθ))
well_coords = [p(k) for k = 0:num_wells-1]
# Plot well coordinates in horizontal plane
scene = scatter(well_coords, markersize=10)

# ## Set up mesh function
# We create a mesh for the BTES system using Gmsh. The mesh is first created in
# 2D with refinement around the wells, and then extruded to create a 3D mesh.
function create_btes_mesh(xw, radius, depth, h_min, h_max, h_z)

    # Set outer radius sufficiently far away from BTES wells
    radius_outer = radius * 5

    # Clear all models and create a new one
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("btes_mesh")

    # Make boundary
    perimeter = 2*π*radius_outer
    n = Int(ceil(perimeter / h_max))
    Δθ = 2*π/n
    for i = 1:n
        x, y = radius_outer * cos(i*Δθ), radius_outer * sin(i*Δθ)
        gmsh.model.geo.addPoint(x, y, 0.0, h_max, i)
    end
    # Line segments along boundary
    for i = 1:n
        gmsh.model.geo.addLine(i, mod(i, n) + 1, i)
    end
    gmsh.model.geo.addCurveLoop(collect(1:n), 1)
    # Top surface
    gmsh.model.geo.addPlaneSurface([1], 1)

    # Add well points
    for (i, x) in enumerate(xw)
        gmsh.model.geo.addPoint(x[1], x[2], 0.0, h_min, n + i)
    end
    nw = length(xw)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, collect(n+1:n+length(xw)), 2, 1)
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h_min, n + nw + 1)

    # Set resolution
    # Distance field
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [n + nw + 1])
    # Threshold field
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", h_min)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", h_max)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 1.1 * radius)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.5 * radius)
    # Background field
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    # Quads instead of triangles
    gmsh.model.geo.mesh.setRecombine(2, 1)

    # Extrude
    num_elements = Int(ceil(depth/h_z))
    gmsh.model.geo.extrude([(2, 1)], 0, 0, depth, [num_elements], [1.0], true)
    gmsh.model.geo.synchronize()

    # Generate mesh
    gmsh.model.mesh.generate(3)
    
    # Read into Jutul mesh format
    mesh = Jutul.mesh_from_gmsh(z_is_depth=true)

    return mesh

end

# ## Create BTES mesh
# The following mesh size parameters will give a mesh with approximately 15k
# cells, and can be adjusted to make the mesh coarser or finer.
# NB: Note that not all mesh sizes will give a valid mesh.
depth = 50.0meter # Depth of the BTES system
# Horizontal minimum/maximum and vertical mesh size
h_min, h_max, h_z = 1.0meter, 25.0meter, 5.0meter
mesh = create_btes_mesh(well_coords, radius, depth, h_min, h_max, h_z)

# ## Create reservoir domain
# Next, we construct a reservoir domain with rock with properties similar to
# that of granite
domain = reservoir_domain(mesh,
    permeability = 1e-6darcy,
    porosity = 0.011,
    rock_density = 2650kilogram/meter^3,
    rock_heat_capacity = 790joule/kilogram/Kelvin,
    rock_thermal_conductivity = 3.0watt/meter/Kelvin,
    component_heat_capacity = 4.278e3joule/kilogram/Kelvin
)

# ## Make model
# We set up a model with a geothermal system and the BTES wells. We model each
# BTES well as a U-tube system, wher the supply and return pipes run parallel
# inside a larger borehole filled with grout. Pipe dimensions and material
# thermal properties have reasonable defaults, and can be adjusted as needed --
# see setup_btes_well for details.
grid = physical_representation(domain)
well_models = []
for (wno, xw) in enumerate(well_coords)
    d = max(norm(xw, 2))
    v = (d > 0) ? xw ./ d : (1.0, 0.0)
    xw = xw .+ (h_min / 2) .* v
    trajectory = [xw[1] xw[2] 0.0; xw[1] xw[2] depth]
    cells = Jutul.find_enclosing_cells(grid, trajectory)
    name = Symbol("B$wno")
    w_sup, w_ret = setup_btes_well(domain, cells, name=name, btes_type=:u1)
    push!(well_models, w_sup, w_ret)
end

model, parameters = setup_reservoir_model(
    domain, :geothermal,
    wells=well_models
)

# Inspect model
plot_reservoir(model)

# ## Initial state and boundary conditions
# We impose fixed boundary conditions at the top surface, enforcing a constant
# pressure of 1 atm and a temperature of 10°C.
state0 = setup_reservoir_state(model,
    Pressure=1atm,
    Temperature=convert_to_si(10.0, :Celsius)
)

geo = tpfv_geometry(mesh)
z_f = geo.boundary_centroids[3, :]
top_f = findall(v -> isapprox(v, minimum(z_f)), z_f)
top = mesh.boundary_faces.neighbors[top_f]

bc = FlowBoundaryCondition.(top, 1atm, convert_to_si(10.0, :Celsius))

# ## Controls
# Each BTES well is modelled using two mutlisegment wells, one for the supply
# and one for the return. The supply well is controlled by a rate target, while
# the return well is controlled by a BHP target. The supply and return wells
# communicate through their bottom cell.

# Rate control for supply side
rate = 0.5litre/second
temperature_charge = convert_to_si(80.0, :Celsius)
temperature_discharge = convert_to_si(10.0, :Celsius)
rate_target = TotalRateTarget(rate)
ctrl_charge = InjectorControl(rate_target, [1.0], density=1000.0, temperature=temperature_charge)
ctrl_discharge = InjectorControl(rate_target, [1.0], density=1000.0, temperature=temperature_discharge)

# BHP control for return side
bhp_target = BottomHolePressureTarget(1atm)
ctrl_prod = ProducerControl(bhp_target)

# Assemble forces
control_charge = Dict()
control_discharge = Dict()
for well in well_models
    if contains(String(well.name), "_supply")
        control_charge[well.name] = ctrl_charge
        control_discharge[well.name] = ctrl_discharge
    else
        control_charge[well.name] = ctrl_prod
        control_discharge[well.name] = ctrl_prod
    end
end

# Set up 6 months of charging, followed by 6 months of discharging
month = year / 12
n_step = 6
dt = fill(month, n_step) # Report every month

forces_charge = setup_reservoir_forces(model, control=control_charge, bc=bc)
forces_discharge = setup_reservoir_forces(model, control=control_discharge, bc=bc)
forces = vcat(fill(forces_charge, n_step), fill(forces_discharge, n_step))
dt = repeat(dt, 2)

n_cycles = 10
forces = repeat(forces, n_cycles)
dt = repeat(dt, n_cycles)

# ## Pack simulation case and simulate
case = JutulCase(model, dt, forces, state0=state0, parameters=parameters)
results = simulate_reservoir(case)

# ## Inspect reservoir states
plot_reservoir(model, results.states)

# ## Inspect well solutions
plot_well_results(results.wells)