# # Embedded fractures on an unchanged reservoir mesh
# <tags: Immiscible, Discretizations, StartToFinish>
# EFM clips prescribed fracture surfaces to matrix cells and adds conservative
# connections. The matrix mesh, cell numbers, and well completions are retained.
using Jutul, JutulDarcy

mesh = reservoir_mesh((12, 12, 4), (120.0, 120.0, 40.0))
matrix = reservoir_domain(mesh; permeability = 0.01si_unit(:darcy), porosity = 0.15)
properties = (aperture = 1e-3, porosity = 0.5)
fractures = AbstractEmbeddedFracture[
    RectangleFracture([60, 60, 20], [1, 0, 0], (100, 35); properties...),
    EllipticalFracture([55, 60, 20], [0, 1, 0], (45, 17); segments = 64, properties...),
    PolygonFracture([[30, 20, 20], [90, 20, 20], [90, 80, 20],
                     [60, 60, 20], [30, 80, 20]]; properties...)
]
embedded = embedded_fracture_domain(matrix, fractures)

# ## Choose whether to unify the models
# Separate matrix/fracture models use the existing DFM cross terms. Setting
# unified=true appends fracture cells and uses explicit transmissibilities.
system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = (1000.0, 800.0))
model = setup_embedded_fracture_model(matrix, embedded, system; unified = false)
single_model = setup_embedded_fracture_model(matrix, embedded, system; unified = true)

# ## Solve flow with pressure boundary conditions
# Boundary cells still refer to the original matrix mesh in both models.
bc = [FlowBoundaryCondition(matrix, 1, 110si_unit(:bar); fractional_flow = [1.0, 0.0]),
      FlowBoundaryCondition(matrix, number_of_cells(matrix), 90si_unit(:bar))]
state0 = setup_reservoir_state(model; Pressure = 100si_unit(:bar), Saturations = [0.0, 1.0])
forces = setup_reservoir_forces(model; bc = bc)
case = JutulCase(model, fill(si_unit(:day), 10), forces; state0 = state0)
result = simulate_reservoir(case; info_level = -1)

# ## Inspect geometry and connections
# The physical fracture mesh can be plotted with Jutul's usual plotting tools.
fracture_mesh = physical_representation(embedded)
println("Matrix cells: ", number_of_cells(matrix))
println("Fracture cells: ", number_of_cells(embedded))
println("Fracture connections: ", number_of_faces(embedded))
println("Intersection segments: ", length(fracture_mesh.geometry.intersections))
