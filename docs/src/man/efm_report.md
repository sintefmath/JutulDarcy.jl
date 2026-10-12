# Embedded fracture implementation results

Initial validation on 2026-10-09; subsequent reviews are recorded below.

The new EFM module supports conductive fractures on an unchanged matrix mesh.
`embedded_fracture_domain` clips rectangles, planar ellipses, and concave
planar polygons and generates fracture cells and transmissibilities. Setup supports both
separate matrix/fracture models with DFM mass and energy cross terms and one
reservoir model with fracture cells appended. Original matrix cell indices and
well completions remain valid.

All **432 regression checks passed**: 158 EFM checks, 203 existing DFM checks,
and 71 mesh and single-phase checks. The two-phase example completed, and the
new methods introduce no dispatch ambiguities. Warped-face clipping and
numerically zero-width patches are handled. Fractures on a matrix face replace
its covered direct connection to avoid duplicate flow paths; the input mesh
and domain remain unchanged.

Aligned fracture geometry and transmissibility stencils match DFM to numerical
precision. The aligned two-phase cube matches DFM to roundoff when both models
use the same storage accounting. Coupled and single EFM simulations agree in
pressure, saturation, and temperature with tightened nonlinear tolerances.

On the nonaligned three-plane cube, DFM states were volume-averaged onto the
original EFM cells. Maximum RMS differences over twenty report steps were:

- **Two-phase, 0.5 pore volumes injected:** relative pressure error was 0.067%
  on the 5×5×3 grid and 0.059% on the 11×11×7 grid. Saturation RMS was 0.0729
  and 0.0833 respectively, so the maximum transient saturation error did not
  decrease under this refinement.
- **Geothermal, 5 pore volumes injected:** temperature RMS decreased from
  1.44 K to 0.637 K on the same grids. Relative pressure RMS decreased from
  0.429% to 0.341%.

These comparisons use matching matrix well completions. EFM retains matrix
cells spanning fractures, while DFM splits them; their nonaligned transport
discretizations therefore differ. Close pressure agreement does not imply
equivalent saturation transport.

Both full active corner-point grids completed preprocessing and flow solves
with **1,200 mixed fracture objects**, using seed 1984:

| Grid | Matrix cells | Fracture cells | Preprocessing | Single solve | Coupled solve |
|:--|--:|--:|--:|--:|--:|
| OLYMPUS | 192,750 | 140,830 | 8.20 s | 4.83 s | 5.36 s |
| NORNE | 44,927 | 93,346 | 7.11 s | 2.28 s | 2.75 s |

These timings cover three one-day single-phase report steps with opposing
110/90 bar pressure boundaries, after compilation, on an Apple M4 with Julia
1.13.1, one Julia thread, one BLAS thread, and HYPRE BoomerAMG. Both
representations produced identical final pressure vectors. Grid loading and model setup
are excluded from the table. NORNE uses its local PERMX values with
PERMY=PERMX and PERMZ=0.1PERMX; this benchmark uses simplified fluid and boundary
conditions rather than the production deck's full setup.

The main remaining performance opportunity is reducing temporary geometry
allocations: preprocessing allocated 13.4 GB cumulatively for OLYMPUS and
10.0 GB for NORNE. These figures are not peak memory usage. The current method
targets conductive fractures and does not model impermeable barriers or
fracture surfaces with holes. Automatic fracture completions are available for
multisegment wells through `extend_wells=true`.

The [EFM guide](embedded_fractures.md) documents the API, limits, and benchmark
reproduction commands.

## Code reuse review (2026-10-12)

EFM and DFM now share unified graph/property assembly, property override
validation, domain copying and completion conversion. Fracture intersection
and well searches share patch-bound construction. EFM uses Jutul's existing
face adjacency and edge–plane intersections alongside its cell bounds, tensor
expansion, cut-cell and DFM utilities.

The common domain builder fixes dropped non-entity metadata in unified EFM and
copies it without aliasing input data. Jutul's edge–plane helper now uses a
scale-relative check, fixing incorrect midpoint intersections on small
Float32 meshes. EFM retains tolerance-aware concave clipping and winding
resolution because the existing convex clipping routines do not cover these
cases. Broader sharing should first establish those geometry contracts in Jutul.

All 16 EFM types and their 64 stored fields are documented, including units,
index spaces, aperture scaling and shared mesh ownership. The guide also
documents domain fields, cut-cell mappings and well completion data.

Validation passed: the complete run had 586 DFM and 640 EFM checks. A fresh
run using the corrected Jutul helper passed 1,528 cut-cell checks and 28
targeted EFM checks, including an end-to-end small Float32 geometry comparison.
Julia registered all added type and boundary API docstrings successfully.

## Internal fracture grids (2026-10-12)

Fracture constructors now accept `resolution=(nu, nv)` or `cell_length=h`.
Preprocessing also accepts global or per-fracture overrides. These create
independent planar grids that can be finer or coarser than the matrix mesh.
Rectangles and ellipses follow their configured axes; polygons use their local
frame. Omitting the options preserves matrix clipping, and `resolution=:matrix`
explicitly restores it. DFM keeps its conforming mesh.

Coarse cells use separate overlap areas and storage weights for each matrix
cell. Disconnected active regions remain separate cells. The implementation
reuses polygon clipping, interface construction and intersection elimination.
The new settings, grid geometry and patch weights are documented in the guide.

All 887 EFM, 586 DFM and 1,534 Jutul cut-cell regression checks passed. New
checks cover independent cell counts, analytical transmissibilities, unequal
overlap weights, preserved area/storage, concave polygons, disconnected regions,
rotation, precision and two-phase/geothermal agreement between formulations.
A shared Jutul adjacency fix also permits meshes containing only isolated cells.

Both full corner-point grids completed three one-day flow steps with 1,000
seeded mixed fractures and `resolution=(4, 8)`. Unified and coupled pressure
vectors agreed within the benchmark's relative tolerance of `1e-5`.

| Grid | Matrix cells | Fracture cells | Preprocessing | Unified solve | Coupled solve |
|:--|--:|--:|--:|--:|--:|
| OLYMPUS | 192,750 | 32,456 | 11.20 s | 16.37 s | 12.81 s |
| NORNE | 44,927 | 31,139 | 8.44 s | 2.96 s | 3.74 s |

These validation timings use the locally cached grids, simplified rock/fluid
properties, opposing 110/90 bar pressure conditions, CSC assembly and HYPRE
AMG. The first OLYMPUS solve includes compilation; a separate regression process
was also running. Preprocessing allocated 16.4 GB and 13.0 GB cumulatively,
respectively, including temporaries; these are not peak memory figures.
Reproduce through `benchmark_embedded_fractures(...; count=1000, resolution=(4, 8))`
in `examples/fractures/efm_cornerpoint_benchmark.jl`.
