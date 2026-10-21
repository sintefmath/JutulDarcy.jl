# 010 — Final retained changes and results (2026-10-21)

## Branches and retained source

All commits and pushes target only `sol-perf` in the original Jutul and
JutulDarcy checkouts. Jutul2/JutulDarcy2 are untouched. The existing local edit
to `examples/validation/5_validation_olympus_1.jl` is excluded from every commit.

Jutul production commits:

- `c7ea1bd3`: fix loop-length minimum batching, respect context worker limits,
  preserve those limits in derived contexts, and use Polyester for KA CPU
  closure loops and AMG CPU row/setup loops.
- `d24bcbcc`: host row batches for scaled CSR products, with dimension and
  beta-zero handling. Accelerator dispatch is unchanged.
- `b648458c`: make the focused worker-limit checks respect Polyester's default
  physical-core cap on SMT hardware; no further production change.

JutulDarcy production commits:

- `aa2fe00e`: specialize and batch the CPU CPR pressure RHS; select existing
  level-scheduled ILU0 for large CPU KA models with default smoother options.
- `c03425a8`: retain serial system ILU below four workers. Automatic parallel
  ILU requires at least 100,000 cells and four default-pool workers. Explicit
  smoother choices/options take precedence. At two workers serial ILU improves
  both fixtures; at four workers Norne is effectively tied (19.222 s serial,
  19.362 s parallel), so keep the parallel choice from four onward.

No new global BLAS setting, worker QoS policy, vector-storage abstraction,
GS schedule/cache, or additional CPR matrix/weight implementation is retained.

## Final controlled scalability

Three one-day report intervals; medians of three fresh simulations after
compilation warm-up; input/simulator construction excluded; initial AMG setup
included. Float64, GMRES, fixed pressure SPAI0, BLAS=1. The system ILU choice
matches the final automatic policy. The two-worker entries were measured with
an explicit serial-ILU override, which the final configuration checks confirm
is now the default at that worker count.

| Workers | Olympus z4, 771,000 cells, s | Scaling | Norne z8, 355,444 cells, s | Scaling |
|---:|---:|---:|---:|---:|
| 1 | 23.022 | 1.00x | 39.232 | 1.00x |
| 2 | 14.984 | 1.54x | 26.278 | 1.49x |
| 4 | 10.489 | 2.19x | 19.362 | 2.03x |
| 6 | 8.995 | 2.56x | 15.773 | 2.49x |
| 8 | 7.490 | 3.07x | 13.706 | 2.86x |
| 10 | 7.143 | 3.22x | 13.024 | 3.01x |
| 14 | 12.731 | 1.81x | 22.656 | 1.73x |

Ten workers is the best measured choice on this M4 Pro (10 performance +
4 efficiency cores). Worker count is not an affinity guarantee. The clean
whole-simulation **4x target was not reached**. The 1,542,000-cell Olympus
experiment gives 3.16x and is additionally affected by memory compression;
it does not establish a larger-size 4x result.

On original grids, the current ten-worker fixed-SPAI0 medians are 0.415 s for
Norne and 0.596 s for Olympus, versus 0.518/0.822 s in the original code:
roughly 20% and 27% lower runtime. These are code-improvement comparisons,
not thread-scaling ratios.

## Choose the pressure smoother explicitly when comparing algorithms

Norne z8 favors GS: 11.899 s at one worker and 5.215 s at ten. It completes
with 17 Newton / 82 linear iterations and no recovered cuts, versus the much
larger fixed-SPAI0 counts in 007. This is 2.28x GS-to-GS scaling, and 2.50x lower
ten-worker runtime than SPAI0. Olympus z4 favors SPAI0 at ten workers:
7.143 s versus 10.175 s with serial GS. The actual serial Olympus GS default
is 19.464 s, making default-to-default scaling about 2.72x. Do not confuse
these comparisons with fixed-algorithm scaling.

Keep the general pressure-smoother default unchanged. The README supplies the
explicit GS command for the short refined Norne fixture. This avoids selecting
one global pressure algorithm from two different convergence behaviors.

## Final GS investigation: leave outside the library

A process-local schedule preserves the natural sweep order on every structural
edge in both directions, including asymmetric/upper-only/lower-only patterns.
All 12 directed-pattern checks pass bit for bit, and the full simulations keep
the GS reference counts and states. Fine-grid schedules have 334 levels on
Olympus (median width 2,099.5) and 461 on Norne (median width 947).

| Ten-worker GS implementation | Olympus z4, s | Norne z8, s |
|---|---:|---:|
| Existing serial sweeps | 10.175 | 5.215 |
| Dependency levels, validate each sweep | 8.810 | 5.080 |
| Batch-size microbenchmark choice | 9.299 | 5.259 |
| Validate schedules at smoother updates | 8.287 | 5.027 |

Batch sizes 8/16/32 provide little isolated benefit; 64/128 regress. The best
combined microbenchmark advantage is below 1%, and the full-run batch variant
regresses. Moving validation to updates helps Olympus somewhat, but still
trails its existing SPAI0 choice. Norne gains only about 0.19 s over serial GS.
Do not add the larger schedule/state/cache implementation on this evidence.
A generic closure-based Polyester wrapper also hit an Apple ARM cfunction
limitation; explicit batch bodies pass. No workaround enters production.

## Verification

- Existing Jutul threading suite plus the focused worker-limit checks: 23 checks
  at one worker and 49 at ten, rerun after the physical-core portability fix.
- Final CPU KA solver policy: six assertions each at two and four workers,
  including cell cutoff and explicit serial/parallel overrides.
- Scaled CSR products: 24 checks at one and ten workers.
- Existing complete Jutul KA context/AD/multimodel suite: one and ten workers.
- Existing complete KA preconditioner suite: ten workers, CPU/JLArrays,
  including native ILU versus serial CSR ILU.
- Existing complete reservoir KA suite: ten workers, including compositional
  wells, hybrid models, reservoir assembly, and the SPE1 adjoint solve.
- Complete final pressure/saturation comparisons across both full worker
  sweeps, original-grid checks, solver variants, the largest fixture, and GS
  prototypes. Counts agree within each fixed pressure algorithm. The largest
  pressure difference is 0.000213 Pa on refined Norne; the largest saturation
  difference is 9.31e-11. Small original/medium comparisons remain below
  1.1e-6 Pa. GS versus SPAI0 is a separate convergence experiment.

Only the small batching/scaled-product regression blocks are added to the
repository tests. The crossover assertions and rejected prototypes remain
local investigation scripts. Raw TOML results, profiles, process logs, and
serialized final arrays remain in ignored `perf_testing/results/`; the
numbered Markdown files preserve the remotely reviewable record.
