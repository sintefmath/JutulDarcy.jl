# Embedded fractures

Embedded fracture models (EFM, also called EDFM) represent conductive planar
fractures without cutting the matrix mesh. The preprocessor clips each surface
to the active matrix cells, generates a lower-dimensional fracture mesh, and
computes matrix–fracture and fracture–fracture transmissibilities. Existing
matrix cell indices, faces, and well completions remain valid.

```julia
using Jutul, JutulDarcy

mesh = reservoir_mesh((20, 20, 5), (200.0, 200.0, 50.0))
matrix = reservoir_domain(mesh; permeability = 0.01si_unit(:darcy), porosity = 0.1)
fractures = AbstractEmbeddedFracture[
    RectangleFracture([100, 100, 25], [1, 0, 0], (150, 40); aperture = 1e-3),
    EllipticalFracture([90, 100, 25], [0, 1, 0], (60, 20); aperture = 1e-3)
]
embedded = embedded_fracture_domain(matrix, fractures)
system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
model = setup_embedded_fracture_model(matrix, embedded, system)
```

The same preprocessed domain supports two formulations:

| `unified` | Reservoir cells | Matrix–fracture coupling |
|:--|:--|:--|
| `false` (default) | Original cells in `:Reservoir`; fracture cells in `:Fractures` | DFM mass and energy cross terms |
| `true` | Original matrix cells followed by fracture cells in `:Reservoir` | Explicit transmissibilities in one connection graph |

```julia
single = setup_embedded_fracture_model(matrix, embedded, system; unified = true)
# The DFM setup entry point also dispatches to EFM for an embedded domain:
coupled = setup_fractured_reservoir_model(matrix, embedded, system)
```

Both forms work with the usual state, forces, wells, and simulation routines.
The unified model's domain stores `domain[:is_fracture, Cells()]` as a
`Vector{Int}`, with `0` for each matrix cell and `1` for each fracture cell.
For the unified formulation, `plot_reservoir(case, result)` displays
matrix and fracture cell values in simulation order. The mesh can also be
converted explicitly with `UnstructuredMesh(physical_representation(reservoir_domain(case)))`.
This copies the original matrix geometry and adds three-dimensional fracture
prisms with their physical aperture, preserving the connection graph as NNCs.
The prisms overlap the matrix; simulation continues to use the EFM domain's
storage and transmissibilities.

The coupled form supports single-phase, immiscible, and geothermal systems. The
single form can use other fluid systems supported by `setup_reservoir_model`.
Wells retain matrix perforations. For an explicit fracture completion in the
single form, first obtain `combined_fracture_domain(matrix, embedded)`, then
pass the appended fracture cell indices to `setup_well` on that domain.

## Use the same fracture objects with EFM and DFM

`setup_fractured_reservoir_model` handles both methods from the same vector:

```julia
efm = setup_fractured_reservoir_model(matrix, fractures, system;
    method = :efm, wells = wells)
dfm = setup_fractured_reservoir_model(matrix, fractures, system;
    method = :dfm, wells = wells)
```

The default `method=:dfm` keeps the cut-cell setup. EFM accepts
`unified=false` or `unified=true`, and DFM supports both options as well. Unified
DFM keeps cut matrix cells first and appends the DFM fracture and retained
intersection cells, using the same flow and heat coefficients as coupled DFM.
Both methods store the integer `:is_fracture` indicator in the unified domain. Any
`AbstractVector` is accepted, including `Vector{Any}`, abstractly typed vectors,
views and vectors with nonstandard indices. Entries are validated as fracture
objects, and fracture IDs always follow the vector's iteration order.
Both methods subtract fracture storage by default through this shared entry
point. Set `subtract_fracture_volume=false` to retain the full matrix storage.
The separate `setup_embedded_fracture_model` entry point remains available.

The DFM path cuts the matrix along fracture planes and polygon boundaries, then
creates the usual `Jutul.EmbeddedMesh` and DFM cross terms. Rectangles, the
polygonal ellipse, and concave polygons keep their finite extent. Auxiliary
cuts refine the matrix around fracture tips without adding fracture cells.
`min_cut_fraction=0.0` is the default to retain these exact boundaries.
`merge_faces=false` preserves the cut-face mappings.

For direct access to domains and mappings:

```julia
dfm_domains = discrete_fracture_domain(matrix, fractures)
cut_matrix = dfm_domains.matrix
dfm_fractures = dfm_domains.fractures
parents = dfm_domains.info[:cell_index]
initial_pressure_dfm = initial_pressure_efm[parents]
```

`info[:cell_index]`, `:face_index`, and `:boundary_face_index` map the cut
entities to their original matrix entities; new faces use parent ID `0`.
`info[:fracture_id]` maps DFM fracture cells to input object IDs, and
`info[:fracture_faces]` gives their cut-matrix faces in fracture-cell order.

Cell properties, face multipliers, fault tags, and depth orientation are
transferred. Explicit transmissibilities and overrides are divided between
split faces by area; new faces use geometric TPFA coefficients. Existing well
completions are distributed over their cut-cell
descendants by volume, preserving the original matrix completion's total well
index. The object-based DFM setup defaults to `extend_wells=false`, matching
EFM's matrix completions; use `extend_wells=true` to enable automatic
fracture completions for multisegment wells in either method and formulation.
EFM intersects well segments with the finite clipped fracture surfaces,
including segments within a fracture plane. Shared patch boundaries do not
duplicate completions, and simultaneous fracture hits share a well node.
Simple wells retain their existing matrix completions. Cell-indexed initial states,
sources and boundary conditions must be
transferred with `parents` when creating a case manually.

Pass preprocessing and per-object property overrides in `fracture_arg`, for
example `fracture_arg=(aperture=1e-3,)`. The direct DFM routines require a
matrix without NNC faces and fractures with interior matrix connections;
boundary-only fractures are rejected. Other face coefficients must be assigned
after cutting. Cut-cell DFM can introduce many
more matrix cells; EFM retains the original mesh.

See `examples/fractures/compare_efm_dfm.jl` for a common water-injection case.

## Fracture collections

Preprocessing groups fractures by concrete type in a `FractureCollection`.
Each member of `collection.groups` is a vector with a concrete element type;
`collection.indices` maps each group to the original fracture IDs. Different
floating-point types form separate groups. Frame construction, property
extraction and clipping run through specialized functions for these vectors,
while collection indexing and fracture-cell ordering retain the input order.

```julia
collection = FractureCollection(fractures)
collection[1] # the first input fracture
collection.groups
collection.indices
```

Both model setup paths build the collection automatically. A prepared
collection can be reused directly. Treat its group vectors as read-only and
rebuild the collection after changing the input list. Collection, geometry and
mesh displays summarize type counts, cells and connections without printing
their coordinate and connectivity arrays.

## Geometry types and precision

### Internal fracture mesh resolution

Each fracture can select an independent internal grid using either
`resolution` (cell counts) or `cell_length` (target spacing in metres):

```julia
rectangle = RectangleFracture([100, 100, 25], [1, 0, 0], (150, 40);
    resolution = (30, 8))
ellipse = EllipticalFracture([90, 100, 25], [0, 1, 0], (60, 20);
    cell_length = 5.0)
embedded = embedded_fracture_domain(matrix, [rectangle, ellipse])
```

`resolution=n` sets `(n, n)` counts. `cell_length=(hu, hv)` specifies different
target lengths along the two axes. Rectangle and ellipse grids follow their
configured axes, including `direction` and `angle`; polygon grids use their
local `FractureFrame` axes. The grid covers the full fracture's bounding box
before clipping to the fracture boundary and active matrix. Target lengths
give `ceil(extent/length)` cells along each axis, with actual spacing no larger
than the target. Ellipse `segments` controls the boundary approximation
independently of the internal mesh.

Preprocessing accepts the same options, either globally or as per-fracture
vectors. Use tuples for axis pairs and vectors for one entry per input object:

```julia
embedded = embedded_fracture_domain(matrix, [rectangle, ellipse];
    resolution = [(20, 6), nothing], cell_length = [nothing, 3.0])
model = setup_fractured_reservoir_model(matrix, [rectangle, ellipse], system;
    method = :efm, fracture_arg = (cell_length = [4.0, 2.0],))
```

Specify one option per fracture. An override replaces that object's mesh
settings; `nothing` retains them. Omitting both options preserves the original
matrix-clipped mesh. Set `resolution=:matrix`, globally or in a per-fracture
vector, to explicitly restore matrix clipping for configured objects.

Internal cells can be finer or coarser than the matrix cells. Coarse fracture
cells exchange with every intersected matrix cell through its own overlap area.
Optional fracture-volume subtraction uses overlap weights, including shared
storage for face-coincident fractures. Clipping can remove cells or split a
grid cell into disconnected active regions; these regions remain separate.
Grid spacing must exceed the geometric tolerance. Cut-cell DFM ignores these
EFM internal grid settings and keeps its conforming mesh.

Effective grids are stored in `physical_representation(embedded).geometry.grids`,
in input fracture order. Each `FractureGrid` records its actual `resolution`,
`spacing`, and local axes; `nothing` indicates matrix clipping. Each patch
stores `matrix_areas` and `matrix_weights` in its `matrix_cells` order.

### Geometry storage

`EmbeddedFractureMesh{T}` stores the original matrix mesh, fracture objects,
tags, and one `EmbeddedFractureGeometry{T}`. Clipped patches, local frames,
intersections, exchange data and the TPFA geometry live in `mesh.geometry`;
the TPFA geometry is accessed with `tpfv_geometry(mesh)`. For example, clipped
surfaces are in `mesh.geometry.patches`.

Preprocessing defaults to the matrix mesh's floating-point type. Set
`embedded_fracture_domain(matrix, fractures; float_type=Float32)` to select it
explicitly. Fracture constructors infer their type from coordinates and also
accept `float_type`. Float32, Float64 and BigFloat geometry are tested. The
default clipping tolerance accounts for coordinate precision; resolving very
small apertures or geometric details on a large coordinate scale requires
sufficient precision in the input coordinates.

The type docstrings below describe every stored field, its units and its index
space. In particular, `FractureFrame.parts` describes the whole input surface,
while `FracturePatch.parts` describes one connected clipped cell. Both use the
frame's local 2D coordinates. `FractureIntersection.cells` and `.faces` use
fracture-local indices. Matrix–fracture exchange endpoints use separate local
matrix and fracture cell indices; unified setup adds the matrix cell count to
fracture indices.

### Domain fields

`EmbeddedFractureGeometry` stores aperture-independent surface geometry.
The EFM fracture `DataDomain` stores physical properties and aperture-scaled
geometry. These are separate so that aperture overrides can reuse the clipped
mesh. Changing a property after domain setup requires recomputing the
coefficients; assigning a new field value alone does not update them.

| Field | Entity | Meaning |
|:--|:--|:--|
| `:aperture` | `Cells()` | Mechanical aperture (m), repeated for each patch belonging to an input object |
| `:fracture_id` | `Cells()` | Input object ID in iteration order; unified matrix cells use `0`. Present for object-based DFM, optional for preprocessed DFM |
| `:cell_normals` | `Cells()` | Unit fracture-plane normal, stored as a 3-by-cell matrix |
| `:volumes` | `Cells()` | Physical bulk volume (m³): patch area times mechanical aperture |
| `:areas` | `Faces()` | Physical fracture connection area (m²): interface length times the mean endpoint aperture |
| `:matrix_face_coverage` | `NoEntity()` | Per-matrix-face covered fraction in `[0, 1]`; used to remove duplicate direct flow and heat paths |
| `:transmissibilities` | `Faces()` | Explicit fracture–fracture flow coefficients (m³), including star–delta elimination |
| `:rock_thermal_conductivities`, `:fluid_thermal_conductivities` | `Faces()` | Explicit heat-conduction coefficients (W/K); the fluid field is scalar per face or has one row per phase |
| `:matrix_cells`, `:connection_cells` | `FractureMatrixConnection()` | Exchange endpoints in original matrix and local fracture cell numbering |
| `:matrix_cell_centroids` | `FractureMatrixConnection()` | Matrix endpoint coordinates (m), 3-by-connection |
| `:matrix_distance`, `:matrix_connection_area` | `FractureMatrixConnection()` | Averaged matrix-to-plane distance (m) and exchange area (m²) |
| `:fracture_matrix_transmissibility` | `FractureMatrixConnection()` | Explicit exchange flow coefficient (m³) |
| `:fracture_matrix_gravity_difference` | `FractureMatrixConnection()` | `g*(z_matrix - z_fracture)` (m²/s²) |
| `:fracture_matrix_thermal_conductivity` | `FractureMatrixConnection()` | Total rock-plus-fluid exchange heat coefficient (W/K) |
| `:is_fracture` | Unified `Cells()` | `Vector{Int}`: `0` for matrix cells, `1` for fracture and retained DFM intersection cells |

Unified cells are ordered `[matrix; fractures]` and faces are ordered
`[matrix; fracture–fracture; matrix–fracture]`. Matrix boundary fields and
non-entity metadata are copied. Other matrix cell fields use the fracture
domain's value when available, otherwise inherit a parent matrix cell value;
fracture net-to-gross defaults to one. Tensor permeability keeps the matrix's
diagonal or six-component symmetric representation. Original domains remain
unchanged.

When `extend_wells=true`, coupled models use the existing DFM well layout.
`well.perforations.fracture` contains fracture-local reservoir cells and
`well.perforations.self_fracture` contains the corresponding well nodes.
Perforation data is stored on `FracturePerforations()` under the original key
with suffix `_frac`, for example `:perforation_centroids_frac` (m),
`:cell_dims_frac` (m), `:well_index_frac` (m³), and
`:thermal_well_index_frac` (W/K). Fracture gravity differences use each
completion's own well node. Unified setup appends these completions to ordinary
`Perforations()`, offsets their reservoir cell indices, and preserves their
computed well indices.

### Reused geometry and utilities

Generic geometry lives in `Jutul.CutCellMeshes`: planar frames and validation,
concave polygon decomposition, winding resolution, tolerance-aware clipping,
cell sections, signed distance integration, segment intersections, spatial
indexing and prism construction. EFM imports these routines; its geometry
code assembles fracture cells, connections and physical coefficients.
`EFM.FractureFrame` is an alias for `Jutul.CutCellMeshes.PlanarFrame` with the
same fields and constructors. The spatial-index types also come from Jutul.

Jutul's cut/glue helpers use the same planar primitives. Their existing
interfaces and default projection axes are preserved. The shared
`plane_basis` selects `convention=:cutting` for those axes or `:axis` for the
public rectangle/ellipse orientation used by EFM. Planar polygon measures use
signed contributions for concave boundaries; warped faces keep Jutul's
face-centre triangle fans. Clipping and normal normalization use scale-aware
checks to preserve small Float32 geometry.

The object-based DFM path calls Jutul's `cut_mesh` and `EmbeddedMesh`. EFM also
reuses Jutul's cell bounds, edge–plane intersections, tensor expansion and
face adjacency. Property override validation, domain copying and unified
graph/property assembly are shared by EFM and DFM inside `EFM`. Well extension
uses DFM's completion and well-index setup; its searches use the same shared
spatial index and clipped-patch bounds as fracture intersections.

## Fracture geometry and properties

`RectangleFracture(center, normal, lengths)` uses two **full side lengths**.
`EllipticalFracture(center, normal, radii)` uses two **semi-axis lengths** and an
inscribed polygon (`segments=64` by default). `EllipsoidFracture` is an alias
for this flattened, planar ellipse, rather than a curved ellipsoid surface.
`direction` sets the first in-plane axis and `angle` rotates the axes in radians.

`PolygonFracture(vertices)` accepts a vector of 3D vertices or a 3-by-n matrix.
Vertices follow the boundary in either direction; concave simple polygons are
supported. Nonplanar, self-intersecting, and zero-area polygons are rejected.
Holes and positive-area overlap between coplanar fracture objects are not
supported. Merge coplanar overlapping surfaces before preprocessing.

All lengths and rock properties use SI units. Object keywords are `aperture`,
`hydraulic_aperture`, `permeability`, and `porosity`. The default permeability is
the cubic-law value `hydraulic_aperture^2/12`. The preprocessor accepts scalar
or per-fracture overrides for these properties. Matrix permeability can be
isotropic, diagonal, or a full symmetric tensor; exchange uses its projection
onto the fracture normal.

## Geometry and transmissibilities

The preprocessor uses a bounding-box tree to locate candidate cells and
intersections. Warped matrix faces use the same face-center triangle fans as
Jutul's unstructured mesh geometry. It supports faults, inactive cells, and
nonconvex planar sections. Fracture connections across matrix cells follow
the matrix's shared-face topology.
Retraced or crossing section boundaries from split, warped face fans are
resolved into positive-winding regions with an exact planar sweep. Disconnected
pieces inside a matrix cell become separate fracture cells.

Matrix exchange uses the exact volume average of absolute distance to the
fracture plane, integrated over signed tetrahedra. This remains finite for
fractures through cell centers. Exchange area is twice the clipped fracture
area for an interior fracture and once per adjacent matrix cell for a fracture
on a cell face. Face-coincident fractures share one fracture cell. Fracture
intersections are split into line intervals and use star–delta elimination,
including multiple fractures sharing the same line. Point contacts do not
produce connections.

Interior fractures retain the matrix's existing transmissibilities. A fracture
coincident with a matrix face replaces the covered fraction of the direct
matrix connection, for both flow and heat conduction, to avoid duplicate paths.
This changes model coefficients on a copied domain; the input mesh and domain
remain unchanged. It is intended for
conductive fractures and does not represent impermeable barriers; use DFM or
a projection-based method for those. Unsplit patches at interior intersections
use averaged distances, so comparison with a conforming DFM mesh is a
discretization comparison. Properties and transmissibilities are fixed by
preprocessing; rebuild the domain after changing geometry or rock properties.
Thermal preprocessing currently accepts one pore-fluid conductivity per cell,
shared by the fluid phases, as in the DFM thermal cross term.

By default, fracture storage is added to matrix storage. Set
`subtract_fracture_volume=true` during model setup to subtract the physical
fracture volume from intersected matrix cells. This copies matrix domain data,
preserves the original mesh, and rejects nonpositive remaining cell volumes.

## Examples and validation

See the [EFM results report](efm_report.md) for the main validation findings
and remaining limitations.

See `examples/fractures/embedded_fractures.jl` for model setup and simulation.
`examples/fractures/efm_cornerpoint_benchmark.jl` provides a seeded benchmark
with 1,200 rectangles, ellipses, and concave polygons on a full corner-point
grid, including preprocessing, both model representations, and flow solves.
Provide local GRDECL, ACTNUM, and rock-property paths; it does not download data.
The benchmark uses scalar pressure assembly and AMG. Set `amg=:jutul` to use
Jutul's native AMG, or `amg=:hypre` for HYPRE.

The tests in `test/efm.jl` compare fracture geometry and connection stencils
with DFM and compare two-phase and geothermal simulation results between the
two representations.

### Recorded validation

Validated on 2026-10-09 with Julia 1.13.1 on an Apple M4, one Julia thread and
one BLAS thread. Both full active corner-point grids used their local PORO/PERM
data and 1,200 mixed fracture objects generated with seed 1984. Timings exclude
grid loading and geometry warmup; simulation methods were already compiled.
For NORNE's PERMX-only property file, the benchmark sets PERMY=PERMX and
PERMZ=0.1PERMX; this is a geometry/flow benchmark, rather than the production
deck's fluid, well, and layer-specific permeability setup.
The flow problem is single-phase, with opposing 110/90 bar pressure boundaries
and three one-day report steps, using scalar CSC assembly and HYPRE BoomerAMG
(`AggNumLevels=0`, `RelaxType=6`) with unscaled BiCGStab.

| Grid | Matrix cells | Fracture cells | Preprocessing | Single flow solve | Coupled flow solve | Linear iterations |
|:--|--:|--:|--:|--:|--:|--:|
| OLYMPUS | 192,750 | 140,830 | 8.20 s | 4.83 s | 5.36 s | 29 |
| NORNE | 44,927 | 93,346 | 7.11 s | 2.28 s | 2.75 s | 35 |

Both representations completed all report steps with identical final pressure
vectors in these runs. Matrix nodes, neighbors, and volumes were checked for
changes. The preprocessor allocated 13.4 GB cumulatively for OLYMPUS and 10.0 GB
for NORNE, including temporary arrays; these are **not peak memory usage**.
NORNE retains its original degenerate-face geometry, whose invalid matrix
half-transmissibilities are handled by JutulDarcy's usual sanitization.

To reproduce with the local test-deck collection:

```julia
include("examples/fractures/efm_cornerpoint_benchmark.jl")
data = "/path/to/test_decks_jutul_darcy/data"
olympus = joinpath(data, "olympus")
benchmark_embedded_fractures(joinpath(olympus, "OLYMPUS/OLYMPUS_GRID.GRDECL");
    property_paths = [joinpath(olympus, "OLYMPUS_1/OLYMPUS_1_PROP_$(key).GRDECL")
        for key in ("PORO", "PERMX", "PERMY", "PERMZ")])
norne = joinpath(data, "norne_reduced/INCLUDE")
benchmark_embedded_fractures(joinpath(norne, "GRID/IRAP_1005.GRDECL");
    actnum_path = joinpath(norne, "GRID/ACTNUM_0704.prop"),
    property_paths = [joinpath(norne, "PETRO/PORO_0704.prop"),
        joinpath(norne, "PETRO/PERM_0704.prop")])
```

DFM geometry and transmissibilities agree to numerical precision on aligned
single, crossed, and three-plane networks. The aligned two-phase cube also
agrees to roundoff when both methods use the same storage accounting.
Coupled and single EFM pressure, saturation, and temperature agree throughout
their transient simulations with tightened nonlinear tolerances.

For the existing nonaligned three-plane cube configuration, DFM matrix states
are volume-averaged onto the original EFM cells before comparison. Maximum
errors over twenty report steps are:

| System / grid | Relative pressure RMS | Transport RMS |
|:--|--:|--:|
| Two-phase, 5×5×3, 0.5 pore volumes injected | 0.067% | Saturation 0.0729 |
| Two-phase, 11×11×7, 0.5 pore volumes injected | 0.059% | Saturation 0.0833 |
| Geothermal, 5×5×3, 5 pore volumes injected | 0.429% | Temperature 1.44 K |
| Geothermal, 11×11×7, 5 pore volumes injected | 0.341% | Temperature 0.637 K |

These are discretization errors between unsplit EFM cells and conforming DFM
cells. Thermal refinement reduces the measured error. The maximum transient
saturation RMS is not monotone under refinement in this case; EFM should not
be treated as an identical transport discretization on cells that straddle a
fracture. The test wells use the same matrix completions in both methods, so
this comparison excludes DFM's automatic fracture well extensions.

```@docs
embedded_fracture_domain
embedded_fracture_mesh
discrete_fracture_domain
discrete_fracture_mesh
EmbeddedFractureGeometry
EmbeddedFractureMesh
EmbeddedReservoirMesh
AbstractEmbeddedFracture
JutulDarcy.EFM.FractureProperties
JutulDarcy.EFM.FractureFrame
Jutul.CutCellMeshes.PlanarFrame
JutulDarcy.EFM.FracturePatch
JutulDarcy.EFM.FractureIntersection
Jutul.CutCellMeshes.BoundingBoxNode
Jutul.CutCellMeshes.CellSpatialIndex
JutulDarcy.EFM.EFMFluidThermalConductivities
JutulDarcy.EFM.FracturePerforationGravityDifference
FractureCollection
JutulDarcy.EFM.FractureMeshOptions
JutulDarcy.EFM.FractureGrid
combined_fracture_domain
setup_embedded_fracture_model
PolygonFracture
RectangleFracture
EllipticalFracture
fracture_polygon
```
