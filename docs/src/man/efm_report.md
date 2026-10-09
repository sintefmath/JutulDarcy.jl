# Embedded fracture implementation results

Validated on 2026-10-09.

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
automatically extend wells into fractures.

The [EFM guide](embedded_fractures.md) documents the API, limits, and benchmark
reproduction commands.
