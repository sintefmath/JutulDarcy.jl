# 007 — Ryzen 7950X continuation (2026-10-21)

## Starting point and objective

Local `sol-perf-2` branches start at Jutul `d24bcbcc` and JutulDarcy
`e9292fbd`. In these Windows checkouts the SINTEF remote is named `upstream`;
the `origin` forks do not contain `sol-perf`. No pushes are authorized.
Earlier numbered results were measured on an Apple M4 Pro; do not combine
their absolute times with this workstation's results.

Target: at least 5x warmed simulation throughput at 15 Julia workers relative
to one worker on this Ryzen 9 7950X (16 cores / 32 hardware threads). Use refined
Olympus and Norne, with the first 20% of report intervals at their original
durations for final validation. Short three-day fixtures are diagnostic only.
BLAS stays single threaded, CPU benchmark processes run sequentially, and
warm-up is separate. Record convergence, iteration counts, phase timings,
allocations, source commits, and complete final pressure/saturation comparisons.
Compare fixed solver settings as well as the application defaults.

Latest user direction: focus on Float64 correctness and execution performance,
including the serial case, aiming for parity with MPI OPM or better. Stop the
mixed-precision and alternative-solver investigation. Keep the existing
nonlinear/timestep controls fixed. Int32 sparse indices still store Float64
coefficients and are a separate memory-traffic optimization. The fivefold goal
has not been achieved, and short screens do not establish OPM parity.

`PERF_FRACTION=0.2` selects the first ceil(20% * report count) intervals and
disables the one-day cap unless `PERF_MAX_DAYS` is explicitly supplied.
The benchmark environment now points at the sibling `.jl` checkouts used here;
MultiComponentFlash is registry-resolved rather than requiring a fourth checkout.

## Investigation queue

1. Reproduce one/15-worker refined baselines and inspect phase profiles.
2. Audit assembly, properties, state updates, convergence reductions, sparse
   products, Krylov operations, CPR/AMG construction and application, ILU,
   wells, allocation/GC, and output for serial work or scheduling costs.
3. Measure small changes individually with the existing interfaces and
   Polyester where appropriate; retain numerical comparisons and rejected probes.
4. Verify complete 20% schedules, intermediate worker counts, correctness
   checks, and commit each supported improvement without pushing.

## Initial short-fixture baseline

Julia 1.13.2 under Ubuntu WSL2, native Linux dependency depot, working-tree
sources from the Windows filesystem. Mesh/input and simulator construction are
outside the timer. Fixed pressure SPAI0 and explicit `:ka_ilu0` system smoother
at both worker counts; this is an algorithm-controlled diagnostic, not a
comparison against the fastest serial application default.

Olympus z2, 385,500 active cells, three one-day intervals:

| Workers | Median simulation, s | Assembly, s | Preconditioner apply, s | Other Krylov, s | Linear setup, s |
|---:|---:|---:|---:|---:|---:|
| 1 | 19.866 | 3.142 | 10.701 | 1.465 | 3.084 |
| 15 | 7.454 | 0.653 | 3.518 | 1.066 | 1.152 |

Scaling is 2.665x. All samples take 22 Newton iterations, 146 linear iterations,
six ministeps, and zero wasted iterations. Complete final pressure and saturation
arrays match bit-for-bit between one and 15 workers for both cases.
Measured 15-worker runs report zero JIT time.

The one-worker call tree identifies sparse ILU row scanning as a large cost:
both triangular solves scan the whole row and test every column despite having
stored diagonal positions and sorted columns. Probe restricted triangular ranges.
Convergence is only about 0.16–0.22 s here, so it cannot by itself close the gap.
Linux exposes sibling CPU pairs 0/1, 2/3, ..., 30/31; compare physical-core
affinity separately if scaling remains limited.

Norne z4 has 177,722 active cells. The corresponding one/15-worker medians
are 23.268 / 8.627 s (2.697x). All runs take 41 total Newton iterations,
143 total linear iterations and seven ministeps, including 15 wasted Newton and
44 wasted linear iterations from recovered cuts. The cost of those cuts is retained.
The measured samples all report zero JIT time.

## First production fix: CPU matrix worker metadata

Jutul `458f2ab9` preserves the CPU context's worker count and Polyester policy
when adapting CSR matrices. Previously transfer unconditionally set one worker
and serial scheduling, which disabled the existing partitioned CSR ILU even
when the caller explicitly requested `smoother_type=:ilu0` at many workers.
Accelerator matrices retain their previous one-worker host metadata.

Validation: direct three-worker context/three-factor ILU check passes inside a
15-worker process. The complete existing Jutul KA context, generic-AD, grouped
multimodel, interpolation, and sparse-product file passes (316 assertions).
The application still automatically selects level-scheduled ILU for large KA
CPU models; benchmark explicit partitioned ILU before changing that policy.

## Existing-solver screen (one warmed sample per configuration)

These are screening samples; repeat finalists and validate full states before
calling a policy change an improvement. Source is Jutul `458f2ab9`.

| Case | Workers | Krylov | System ILU | Pressure smoother | Seconds | Newton / linear |
|---|---:|---|---|---|---:|---:|
| Olympus z2 | 1 | GMRES | serial CSR | default GS | 8.625 | 14 / 70 |
| Olympus z2 | 1 | BiCGStab | serial CSR | default GS | 8.464 | 14 / 42 |
| Olympus z2 | 15 | GMRES | partitioned CSR | SPAI0 | 7.441 | 22 / 151 |
| Olympus z2 | 15 | GMRES | level scheduled | SPAI0 | 7.401 | 22 / 146 |
| Olympus z2 | 15 | BiCGStab | partitioned CSR | SPAI0 | 6.333 | 16 / 84 |
| Olympus z2 | 15 | BiCGStab | level scheduled | SPAI0 | 5.526 | 15 / 69 |
| Norne z4 | 1 | GMRES | serial CSR | default GS | 20.037 | 41 / 163 |
| Norne z4 | 1 | BiCGStab | serial CSR | default GS | 26.336 | 43 / 175 |
| Norne z4 | 15 | GMRES | partitioned CSR | SPAI0 | 8.880 | 38 / 200 |
| Norne z4 | 15 | GMRES | level scheduled | SPAI0 | 8.469 | 41 / 143 |

Partitioned ILU is now available through its existing interface but has not yet
shown an end-to-end advantage on these short cases. Do not replace the default
with it on this evidence. The true serial default is much faster than fixed
SPAI0/level-ILU on Olympus; the earlier 2.665x ratio overstates scaling relative
to the best serial setting. Pressure-smoother strength and Krylov choice change
Newton counts, so investigate these together rather than extrapolating kernel
speedups. Norne's recovered cuts remain included (see raw TOML).

Initial long workload: original-duration 20% schedules on Olympus z4 (771,000
cells) and Norne z8 (~355,444 cells), one warmed screening sample at 15 workers,
with `PERF_CPUSET=0-28:2` to use distinct physical cores. Olympus uses BiCGStab;
Norne uses GMRES. Both retain current level ILU/SPAI0. This establishes the
longer-workload phase distribution before selecting further code changes.

## Probe corrections

An initial synthetic seven-point matrix retained wraparound edges, producing
192,000 dependency levels of width one. Its ILU times are invalid as a
parallel-grid model. The corrected 80x60x40 grid removes those edges and has
178 levels. One-worker diagonal-bounded traversal improves scalar / 2x2 / 3x3
apply from 3.637 / 12.126 / 19.251 ms to 3.116 / 10.270 / 17.378 ms, with
bitwise-equal complete solutions and zero allocation. At 15 workers, however,
the candidate takes 0.862 / 1.681 / 3.359 ms versus 0.795 / 1.532 / 3.696 ms
for the original, and allocates 22,528 bytes per application. Reject this version:
it regresses the scalar and 2x2 cases and is not a consistent threading win.

The first reduction prototype captured a reassigned accumulator in a closure,
creating allocations and an unsafe reduction. Moving volume reduction to a
separate function removes that capture. An AD-array multiworker comparison
failed because the synthetic dual used a non-Jutul tag, so `Jutul.value`
correctly retained it and the reduction changed type mid-loop. Use Jutul's
existing tagged-AD constructor in the fixture and keep this out of production
until the corrected multiworker comparison passes.

The corrected tagged-AD probe now passes all cases at one and 15 workers.
At 385,500 cells / 15 workers, the two-phase tagged case improves from
3.280 to 0.202 ms and the three-phase tagged case from 7.327 to 0.322 ms.
Temporary allocation falls from 6.17 / 9.25 MB to 240 / 320 bytes per call.
The serial generic black-oil reduction remains faster on plain Float64 arrays,
so retain that path for one worker. `cpu_convergence_candidate.jl` contains
the isolated helper and additional tagged-AD / NaN / shape checks; it is not
included by package source. Validate complete simulations before integrating it.

## Long schedule screening and pressure smoothers

Olympus z4's selected 20% prefix is four original annual intervals:
366 + 365 + 365 + 365 days. It is much more expensive than the three-day
fixture. At about nine minutes into its warm-up, a single documented Julia
SIGUSR1 profile peek confirmed active preconditioner updates, BiCGStab, ILU,
and assembly, rather than compilation or an idle process. Its call stack points
to the warm-up call. The 1-second peek overflowed its default sample buffer;
use it qualitatively only. Use fresh uninstrumented repeats for final claims.
Process CPU use is roughly 10–12 cores, RSS about 6 GB, with no swap use.
Future long warm-ups log physical progress; measured runs remain quiet.
Stopped this diagnostic while still warming after more than 14 minutes: it
produced no timed sample. The following Norne run did not start. Preserve the
peek and log as evidence of active work, not as a scaling measurement.

The user directed this continuation to focus on KA because HYPRE threads poorly.
Stopped the queued standard-CPU screen; its completed one-worker samples are
retained for reference only. No HYPRE policy or production change is proposed.
The active KA screen compares existing `:gauss_seidel`, `:dilu`, and `:ilu0`
pressure smoothers with both GMRES and BiCGStab, level-scheduled system ILU,
15 physical cores, one BLAS worker, and the same short refined fixtures.
Labelled result files prevent overwriting the earlier solver screen.

The isolated fused host convergence helper passes 27 checks at one worker and
27 at 15, including tagged AD, black-oil and general multiphase equations,
NaNs, dimensions, and exclusion of external AD. It is still outside package
source pending full-simulation before/after validation.

## KA pressure smoother screen

One warmed short-fixture sample per configuration, 15 physical cores, level-
scheduled system ILU, Float64. These are screening results, not final scaling
claims. No defaults have changed.

| Case | Krylov | Pressure smoother | Seconds | Total Newton / linear | Wasted Newton / linear |
|---|---|---|---:|---:|---:|
| Olympus z2 | BiCGStab | GS | 4.919 | 14 / 42 | 0 / 0 |
| Olympus z2 | BiCGStab | ILU0 | 41.425 | 20 / 389 | 0 / 0 |
| Olympus z2 | GMRES | GS | 4.868 | 14 / 70 | 0 / 0 |
| Olympus z2 | GMRES | ILU0 | 11.749 | 15 / 172 | 0 / 0 |
| Norne z4 | BiCGStab | GS | 11.814 | 46 / 124 | 15 / 28 |
| Norne z4 | BiCGStab | ILU0 | 15.599 | 48 / 130 | 15 / 29 |
| Norne z4 | GMRES | GS | 9.635 | 41 / 163 | 15 / 70 |
| Norne z4 | GMRES | ILU0 | 8.419 | 29 / 97 | 0 / 0 |
| Norne z4 | GMRES | DILU | 58.655 | 57 / 1328 | 30 / 973 |

Olympus DILU exceeds the 120-second per-run screening budget with both Krylov
methods, with repeated cuts and poor linear convergence. Norne BiCGStab/DILU
also exceeds that budget. The first exploratory
Olympus DILU warm-up was stopped manually before adding this bound. A single
SIGUSR1 peek during that warm-up confirmed active solver work. No timing from
the instrumented warm-up is used above.

GS is best on Olympus. Norne does not share that outcome: ILU0/GMRES avoids
recovered cuts but has a higher per-iteration cost, finishing near the earlier
SPAI0/GMRES time. Test the original-duration schedules and AMG update policies
before changing application defaults. The shared comparison driver now uses
the existing iteration hook to bound candidates, and forwards named options to
the existing solver API. An unbounded driver remains available for final runs.

The broader preconditioner audit found two independent native CPU correctness
issues: ILU `apply!` ignores non-default alpha/beta and loses the old output,
while CSC SPAI0 fails to clear its normalization buffer on numeric updates.
Focused regression tests are prepared. These are outside the current KA hot
path and do not explain the smoother timings above. Keep their fix separate
from performance claims.

## Retained host convergence reduction

Fused reductions are now integrated into the existing multiphase and black-oil
functions for large strided Float64 host arrays at multiple workers. Other
arrays, external/nested AD, non-Float64 values and serial runs retain the
generic code. No exported API is added. The prototype was removed, and its
tests now exercise both the private reduction and the production entry points
against independent serial references: 46 checks at one worker, 54 at 15.

Three fresh warmed simulations with the original controlled SPAI0/GMRES/level-
ILU configuration and unrestricted affinity give 6.945 s for Olympus z2 and
7.937 s for Norne z4, versus 7.454 / 8.627 s previously. All six measurements
report zero compilation time. Complete final arrays remain bit-for-bit equal
to the original 15-worker baseline, with unchanged total/wasted counts.
These are observed improvements of about 6.8% / 8.0%, not fivefold scaling.
Frozen pressure/system fixtures are saved outside timed simulation work for
focused preconditioner replays; their IO is excluded from the timing.

Jutul `3baad4b2` fixes the two native CPU preconditioner bugs independently.
New tests reproduce nine failures before the fix and pass all 24 checks after
it at both one and 15 workers. Existing scalar/block CSR ILU tests also pass
all 24 checks at both worker counts. The alpha=1, beta=0 fast path is retained.

Incorporated Jutul's portable batching-test commit `b648458c` and the other
session's final M4 report (`2d5ee982`, local `332f1a72`). Its GS prototypes
preserve states but provide small end-to-end gains. The M4 system-ILU crossover
commit `c03425a8` is deferred here: after the worker-metadata fix, native `:ilu0`
at two workers is partitioned rather than serial. Validate that policy on Ryzen
before importing a threshold whose stated serial behavior no longer applies.

## Original-duration 20% screening

Olympus z4 completes four original annual reports (1,461 days) in 857.859 s at
15 physical cores with GS pressure AMG, level-scheduled system ILU and GMRES.
This is one diagnostic sample with progress output, after a separate three-day
warm-up. It reports 3.034 s compilation (0.128 s recompilation) and is not
a final scaling measurement. Final pressure/saturation arrays are retained.

The counters are totals including failed attempts: 1,517 Newton / 5,829 linear
iterations, of which 840 / 2,319 are wasted, over 178 attempted ministeps.
Earlier descriptions that called the totals accepted work were incorrect;
the raw TOML counters and timings are unchanged. Failed work accounts for
55% of the Newton total. At the user's direction, leave timestepping and nonlinear
controls unchanged. Investigate the scaling of each operation and report total
and wasted work separately when solver changes alter iteration counts.

Preconditioner application takes 456.22 s, linear setup/update 128.64 s, and
Krylov work 111.82 s. Together these are about 81% of wall time. The warmed
pressure hierarchy has 13 levels, about 1.44 million total rows, and an operator
complexity around 2.7. Profile/replay pressure smoothing and coarsening as well
as the system smoother. Convergence checks take 10.10 s after the fused fix.

Norne z8's original 20% prefix completes 50 reports, spanning roughly 1.57 years,
in 1,171.912 s with the same settings. It reports 2.254 s compilation (0.127 s
recompilation), 2,834 total Newton / 8,433 total linear iterations, of which
1,200 / 3,160 are wasted, over 363 attempted ministeps. Preconditioner application
takes 489.21 s, setup/update 207.64 s, Krylov work 110.40 s, assembly 179.12 s,
and secondary properties 83.86 s. Both diagnostics use original report durations;
the capped warm-up is recorded separately. Neither has a matched serial 20%
reference yet, so neither establishes a scaling ratio.

New process-local GS tests enforce ordering on every structural edge in either
direction. A forward level order can be reversed for backward GS because both
respect every natural-index adjacency. The row operation remains the existing
GS sweep. All 522 Float64 checks pass at one and two workers, including directed
graphs, damping, cache reuse after numeric changes, and pattern invalidation.
Those tiny correctness checks ran on the physical core excluded from the Norne
diagnostic; final timing repeats will be quiet and run alone. Extended Float32/
Float64 coverage subsequently passes 1,074 checks at both two and 15 workers.
The first pressure replay nevertheless takes 16.83 ms versus 16.36 ms for serial
GS, with allocations rising from 2,896 to 369,008 bytes per application. Do not
retain that implementation. A stable function boundary and fixed global row
ownership are being tested in isolation, including directed graphs with wide
levels and optional per-worker Linux affinity. No schedule/cache implementation
has entered Jutul. The WSL guest exposes one L3 cache domain, so guest CPU numbering
does not establish the host's CCD mapping.

## Further preconditioner screens

At 15 physical cores on Olympus z4, three one-day reports with GS pressure AMG
and GMRES give the following single warmed samples. They screen candidates;
they do not replace the original-duration measurements.

| System smoother | Seconds | Total Newton / linear | Wasted Newton / linear |
|---|---:|---:|---:|
| KA ILU0 | 19.836 | 25 / 157 | 0 / 0 |
| KA DILU | 19.796 | 25 / 157 | 0 / 0 |
| Native partitioned ILU0 | 20.184 | 24 / 175 | 0 / 0 |
| KA SPAI0 | 26.821 | 29 / 291 | 0 / 0 |

DILU's behavior on the system matrix differs from its poor pressure-AMG result.
It is viable here, but does not provide a material measured gain. Norne z8 gives
10.544 s / 17 Newton / 82 linear with KA ILU0 and 9.918 s / 17 / 83 with KA DILU,
both without failed work. Native partitioned ILU takes 27.891 s / 48 / 289
(15 / 78 wasted), and KA SPAI0 takes 95.869 s / 91 / 1,603 (45 / 922 wasted).
These are single warmed samples; final arrays are saved for comparisons.
Leave all timestep/nonlinear controls fixed while evaluating preconditioners.

The CPR initial-copy elimination probe preserves results bit-for-bit but gives
no consistent benefit: the default pre/post 1/0 configuration changes 20.575 to
20.521 ms on Olympus z2 and 13.967 to 14.626 ms on Norne z4. Reject that change.
The corrected, allocation-free ILU range probe also remains inconsistent at 15
workers (scalar 0.732 to 0.830 ms; block-2 1.095 to 1.003; block-3 2.726 to 2.718).
Its serial improvements alone do not justify changing the general path.

The revised isolated GS scheduler passes 2,853 Float32/Float64 checks at both
one and 15 workers. Coverage includes fixed/rank-based ownership, both sweep
directions, correction sweeps, damping, asymmetric edges, and wide levels that
exercise all 15 workers. Pressure-matrix performance replays are in progress;
correctness alone is not a reason to retain the prototype.

## Fixed-input operation scaling

`probe_terms.jl` replays the same frozen arrays at 1, 2, 4, 8 and 15 workers,
each pinned to a distinct physical CPU. BLAS stays at one thread. Measurements
are medians of nine batches after warming each operation, with compilation
recorded separately. These are operation timings, not end-to-end speedups.
The first Olympus serial run preceded addition of the two-sweep SPAI case;
`term_scaling_olympus_t1_v2` repeats it with the final schema. Other runs use
that final schema. Numerical library sources remain `3baad4b2` / `2ffcf2b7`.

| Term, Olympus z2 | 1 worker, ms | 15 workers, ms | Speedup |
|---|---:|---:|---:|
| Block CSR product | 4.778 | 1.766 | 2.70x |
| System ILU application | 35.713 | 5.039 | 7.09x |
| System ILU numeric update | 54.228 | 11.894 | 4.56x |
| Full CPR, one SPAI pressure sweep | 63.780 | 20.426 | 3.12x |
| Pressure AMG, one SPAI sweep | 19.456 | 7.396 | 2.63x |
| Pressure AMG, two SPAI sweeps | 25.272 | 9.993 | 2.53x |
| Pressure AMG numeric update, SPAI | 43.260 | 18.430 | 2.35x |
| Pressure AMG, GS | 22.060 | 17.196 | 1.28x |
| Krylov vector norm | 0.287 | 0.289 | 1.00x |

| Term, Norne z4 | 1 worker, ms | 15 workers, ms | Speedup |
|---|---:|---:|---:|
| Block CSR product | 6.947 | 1.641 | 4.23x |
| System ILU application | 24.810 | 4.799 | 5.17x |
| System ILU numeric update | 37.680 | 13.009 | 2.90x |
| Full CPR, one SPAI pressure sweep | 46.872 | 14.664 | 3.20x |
| Pressure AMG, one SPAI sweep | 10.435 | 2.428 | 4.30x |
| Pressure AMG, two SPAI sweeps | 13.109 | 3.105 | 4.22x |
| Pressure AMG numeric update, SPAI | 22.674 | 8.049 | 2.82x |
| Pressure AMG, GS | 10.265 | 9.329 | 1.10x |
| Krylov vector norm | 0.198 | 0.198 | 1.00x |

AMG residuals/restriction/prolongation often scale 6–16x in isolation, while
serial GS sweeps do not scale. Isolated per-level measurements repeatedly reuse
one level's arrays; their sum is not a prediction of whole-cycle time. Larger
working sets and traffic across the hierarchy are a hypothesis for the remaining
gap. Measure full well-reduced operators and individual numeric-update stages,
and test existing Int32/mixed-precision options separately. No timestep or
nonlinear-control changes are part of this investigation.

## Matching-work full-simulation follow-up

Three warmed Olympus z4 samples (three one-day reports) give 32.740 s with
native serial ILU/GS and 20.056 s with 15-worker KA level ILU/GS. Both use 25
Newton / 157 linear iterations without wasted work; this comparison is about
1.63x. Native serial ILU with two/four SPAI sweeps takes 34.577 / 36.527 s;
the corresponding 15-worker level-ILU times are 19.156 / 18.692 s. The first
two-sweep parallel sample contains 0.668 s compilation following a diagnostic
SIGUSR1 peek during its warm-up; repeat that configuration independently before
using it as a final reference. The GS and four-sweep samples report zero JIT.

For matching-work Olympus GS, median phase times show the remaining limits:

| Phase | Serial, s | 15 workers, s | Speedup |
|---|---:|---:|---:|
| Assembly / linear-system update | 7.055 | 1.530 | 4.61x |
| Secondary properties | 1.328 | 0.487 | 2.73x |
| Equation work outside fused assembly | 0.146 | 0.113 | 1.29x |
| Convergence | 0.383 | 0.143 | 2.68x |
| State updates | 0.407 | 0.263 | 1.55x |
| Krylov excluding preconditioning | 3.620 | 2.441 | 1.48x |
| Linear setup | 5.037 | 2.627 | 1.92x |
| Preconditioner application | 13.889 | 11.667 | 1.19x |

Norne z8's 15-worker medians are GS 9.965 s (17 / 82 Newton/linear), two SPAI
sweeps 15.471 s (28 / 129), and four sweeps 38.036 s (95 / 299, including
60 / 198 wasted). Its pressure-only replay winner therefore does not establish
a full-simulation improvement. Keep separate model-specific results.

The secondary-property path already processes all variables within each worker's
cell interval. Those intervals can be tens of thousands of cells. A possible
future probe is smaller cache tiles inside each worker interval, retaining the
configured worker limit and variable dependency order; simply increasing the
number of dispatched batches would lose the context's worker cap. No source
change has been made for this hypothesis.

## Recycled ILU symbolic state correction

Jutul `c1379eaf` fixes a demonstrated ILU0/DILU setup bug. AMG's
`owned_on_backend_reusing` can replace the contents of retained CSR index arrays.
The old smoother-compatibility check accepted those arrays by identity, reusing
stale diagonal positions and dependency levels. A same-size/same-nonzero-count
graph replacement reproduces 32 failures out of 40 checks; one reused solve
returns `[0.25, -2, 48.75, 1]` instead of approximately
`[0.518, 0.768, 1.071, 1]`.

Setup now checks the saved index snapshots, with allocation-free host comparisons.
The ordinary numeric-update shortcut is unchanged. The regression passes all 40
checks at both one and 15 workers, alongside 24 native preconditioner checks,
24 scalar/block CSR precision checks and 48 isolated hybrid-smoother checks at
each count. Broader KA regressions pass 1,271 assertions across 51 summaries,
including reservoir adjoint checks. Earlier pressure ILU/DILU full-
simulation screens predate this correction and need rechecking before treating
their results as definitive rejections. This correctness fix is not a claimed
speedup of the GS/SPAI baseline paths.

The hybrid experiment uses existing SPAI0 kernels above a row-count threshold
and existing GS kernels below it, through the existing smoother interface.
It adds no application default or public option. Its setup, numeric reset and
memory-rebuild checks pass at one and 15 workers. Pressure-only tests at 15
workers do not establish an advantage over SPAI0: Olympus SPAI0 takes 7.318 ms,
hybrid 32k takes 7.267 ms; Norne takes 2.002 versus 2.591 ms. Full-simulation
performance remains unvalidated. The separate full-operation probe measures
measure saved-hierarchy application, a native serial ILU reference, Schur
couplings, CPR transfers, Galerkin updates and a fixed-system GMRES solve.

Stopped the last Norne serial four-sweep SPAI screening after it consumed several
minutes without completing the comparison. Completed GS/two-sweep samples and
their full arrays are retained. No four-sweep serial timing is accepted for
Norne. Remote branches checked again at 19:36 UTC are unchanged.

## Fixed-input locality and existing operation benchmarks

The property-tiling probe uses the same initial state and verifies every
secondary array bitwise. Olympus z4's original update takes 51.656 ms at one
worker and 19.313 ms at 15; 2,048-cell tiles take 12.438 ms at 15. Norne z8
takes 85.933 / 27.393 ms originally and 20.243 ms with tiles at 15. However,
the generic tiled loop allocates 5.286 / 9.297 MB per call, respectively,
compared with 0.215 / 0.810 MB originally. Do not adopt that implementation.

At the user's suggestion, inspected the existing `Jutul.Benchmarking` routines
for individual secondary variables, assembly and linear operators. A thin
`probe_builtin.jl` wrapper records their timing trees without integrating the
case. Its host per-variable timings differ from production's outer worker
batching; report that distinction. These measurements are queued at one and
15 workers on Olympus z4 and Norne z8.

The full-operation replay initially compared a fresh factor of the saved
assembly matrix against a saved factor from the last linear solve. Those may
represent different Newton iterates: the converged Jacobian can be reassembled
without another linear solve. The comparison now refactors both native and KA
ILU from exactly the same matrix and passes. Olympus native serial apply takes
5.294 ms, while KA takes 33.901 ms at one worker and 5.022 ms at 15. Thus the
apparent 6.75x KA-to-KA application speedup is only 1.05x relative to the fast
native serial implementation. Whole fixed-system GMRES using the saved CPR
takes 1.393 s / 0.529 s; this uses the saved SPAI0 hierarchy, not a GS baseline.
All these measured operations report zero JIT. Index-width and mixed-precision
screens use existing interfaces with fixed nonlinear/timestep controls.

## Dispatch and narrower-index follow-up

Jutul `c1379eaf` / Darcy numerical source `2ffcf2b7`, BLAS one worker,
15 individually pinned default workers. Two fresh full-short samples per
median, all measured JIT zero:

| Case / pressure smoother | Int64 indices, s | Int32 indices, s | Newton / linear |
|---|---:|---:|---:|
| Olympus z4 / GS | 20.056 | 18.332 | 25 / 157 |
| Olympus z4 / four SPAI sweeps | 18.692 | 15.635 | 23 / 135 |
| Norne z8 / GS | 9.965 | 9.506 | 17 / 82 |

Complete GS state comparisons against native serial references give pressure
errors below 5.4e-6 Pa and saturation errors below 4.0e-12. Olympus Int64 serial
and parallel GS are bitwise identical. These improvements do not establish 5x
scaling. The corrected Norne Int32 pressure-ILU screen takes 72.136 s with
142 / 445 Newton/linear iterations; exclude it from remaining precision screens.

The user's launch-overhead observation is material. System ILU application has
548 dependency-level dispatches on Olympus z2 and 586 on Norne z4. Combined
level widths have medians 1522.5 / 649 rows and maxima 2634 / 1114. A store-only
loop with the same dispatch shapes takes 1.476 / 1.156 ms at 15 workers versus
0.036 / 0.017 ms at one. This includes stores and host loop overhead; it is
not a pure launch-latency measurement.

| Fixed saved factors, 15 workers | Olympus z2, ms | Norne z4, ms |
|---|---:|---:|
| Existing ILU application | 4.224 | 4.522 |
| Minimum batch 8 | 4.002 | 3.910 |
| Minimum batch 128 | 6.540 | 7.550 |
| Minimum batch 512 | 12.945 | 22.480 |

All batch variants are allocation-free and bitwise equal to the existing solve.
Repeat favorable choices before changing defaults: minimum batch 1 also differs
from the existing function's first timing, so a single ordered sweep is not
sufficient. The natural-order serial traversal of the same unsplit factor
storage takes 9.784 / 13.517 ms, compared with 33.577 / 24.991 ms through serial
dependency levels. Its relative differences are below 1e-18 / zero. Traversal
locality is therefore a separate cost from thread dispatch. A further probe
uses existing native split-CSR construction and triangular row kernels with the
same dependency schedules, retaining the saved factors rather than refactoring
a potentially different Newton iterate.

The specialized tiled property experiment uses a tuple and `foreach` to keep
variable types and symbols concrete. On Olympus z4 at 15 workers, the original
update is 19.698 ms / 214,832 B, while 512-cell tiles take 11.428 ms / 2,032 B.
Every secondary array is bitwise equal. Serial variants measure zero allocation.
This resolves the generic loop's multi-megabyte allocation growth. Norne and
full-simulation validation are still pending; no library property-loop change
has been adopted. Existing Jutul operation benchmarks are queued next.

The same typed experiment on Norne z8 takes 18.136 ms with 512-cell tiles,
versus 27.255 ms originally, with allocation reduced from 810,192 to 3,776 B.
Every secondary array is bitwise equal at one and 15 workers. A plain named-
tuple `foreach` still allocates megabytes when tiled; explicit `Val` names
are needed for this prototype. The per-call definition-construction cost and
full-simulation behavior remain to be validated. No library source is changed
by these process-local experiments.

## Existing operation benchmarks: Olympus z4

The initial-state built-in benchmark uses Jutul's own routines. Each variable
has 50 repeated evaluations; assembly has 30 and linear products have 100.
Both columns below use its threaded per-cell host path, with one/15 workers:

| Secondary variable | One worker, ms | 15 workers, ms | Speedup |
|---|---:|---:|---:|
| Fluid volume | 1.993 | 0.863 | 2.31x |
| Phase densities | 8.827 | 0.920 | 9.59x |
| Total masses | 4.610 | 2.632 | 1.75x |
| Relative permeabilities | 14.209 | 2.910 | 4.88x |
| Phase viscosities | 15.424 | 1.206 | 12.79x |
| Phase mobilities | 4.901 | 2.336 | 2.10x |
| Phase mass mobilities | 4.118 | 1.983 | 2.08x |

The separate serial whole-variable path has different overhead; e.g. fluid
volume takes 1.036 ms there. Do not silently use that column as the one-worker
threaded baseline. Production still batches all variables within a worker's
cell interval. A further existing-interface probe will compare per-variable
range batches with per-cell dispatch, checking all arrays before timing.

Assembly's linear-system update takes 246.126 / 54.704 ms, and its equation
update excluding secondary variables takes 4.834 / 3.646 ms. Combined this is
4.30x. The full well-reduced operator takes 14.780 / 5.815 ms (2.54x). Repeated
CPR update takes 219.799 / 96.456 ms and application 137.630 / 64.662 ms. These
use KA system ILU even at one worker; the earlier native serial reference is
faster. They measure the fixed initial state and repeated coefficients, not a
whole simulation. The existing `reduce_memory=false` assembly interface is
queued separately to compare compute/storage tradeoffs without new kernels.

The first built-in wrapper completed its numerical work but failed on cell-
count metadata from the stripped evaluation model. It now reads that count
from the full case model before timing and checkpoints each finished stage.
The completed Olympus results above come from the corrected rerun.

## Audit status

| Area | Evidence / next check |
|---|---|
| Fused assembly and secondary properties | Existing Polyester paths; Olympus assembly scales 4.8x. Retain state comparisons. |
| Linear algebra and Krylov | Main products already use host batches. Serial BLAS/vector operations scale weakly. Other session rejected a vector wrapper due to allocation growth. |
| ILU factorization and application | Worker metadata fixed; narrowed row traversal rejected after multiworker regression. |
| Pressure AMG | Host row batches already present; setup and coarse operations remain partly serial. Separate construction from repeated numeric updates in replay probes. |
| CPR weights and pressure matrix | Remaining KA CPU launches; measure host-loop dispatch without changing GPU dispatch. |
| Convergence | Integrated fused reductions pass one/15-worker checks and complete-state comparisons. |
| State updates and pressure increments | Existing strided broadcasts; measure only if profiles justify changes. |
| Wells and Schur complement | Existing per-well blocks and host/device policies; retain complete reservoir KA regression coverage. |
| Allocation, GC, output | Retain GC/allocation counters and saved-state IO in simulation timing; fixture serialization is optional and outside timed runs. |

## Concurrent work

Checked remote `sol-perf` at the user's request. Jutul is unchanged; JutulDarcy
advanced to `bc8be059`, adding M4 larger-grid results and a rejected dense-vector
wrapper probe. Its best reported scaling is about 3.22x on Olympus z4 at ten
workers. Incorporated that commit after the baseline batch. No remote writes.

Subsequently incorporated `8c99c155`, `17658825`, and `ab95c4f1` (documentation).
The latter records no useful full-run gain from CPR host row loops or narrowed
ILU scans, a major threaded-BLAS regression, and a planned correct dependency-
scheduled GS prototype including asymmetric edges. Avoid duplicating that GS
implementation; test resulting commits on Ryzen when available. This continuation
focuses meanwhile on host convergence reductions, KA smoother selection,
original-duration schedules, and physical-core scheduling.

## Float64 execution follow-up

All subsequent simulation and kernel work uses Float64. Historical mixed-
precision screens are evidence only; the queued mixed-convergence experiment
was disabled when the user clarified the scope. The installed local reference
is OPM Flow 2026.04 with MPI. A deck-export helper reuses OPM's own parser and
mirrors `refine_layers.jl`; geometry, property, completion and schedule checks
must pass before its timings can be compared. Final comparison uses the first
20% of the original schedules, recording both setup and simulation time and
iteration counts. Different nonlinear work between programs must remain visible.

Typed property construction per call retains the low allocation counts.
Separately batching every whole variable, or every cell, did not outperform
cache tiling on these large grids: Olympus 20.478 / 19.815 ms and Norne
29.426 / 34.613 ms for the two whole-variable modes, compared with about
11–12 / 18 ms for typed tiles. Repeated isolated-variable benchmarks keep a
smaller working set hot and do not predict the whole update by summation.

Three fresh Float64/Int32 short simulations with 1,000-cell tiles completed at
15 workers. Olympus z4 median is 18.713 s (25 Newton / 157 linear iterations),
Norne z8 9.753 s (17 / 82). Every final pressure and saturation is bitwise
identical to the corresponding untiled Float64/Int32 reference. Secondary
update time falls to about 0.334 / 0.339 s, but total runtime is not consistently
better than the earlier 18.332 / 9.506 s references. No production tiling change
is justified by that unpaired comparison. `compare_execution.jl` now alternates
baseline, tiles, raw GS, and both in one process with fresh simulators, reverses
every second round, and saves complete-state differences and phase timings.

The raw-array dependency-scheduled GS prototype passed 2,853 exact-result tests
at both one and 15 workers, including directed/asymmetric graphs and wide
levels. Its first large pressure replay is:

| Float64 pressure cycle at 15 workers | Original GS, ms | Raw arrays, ms | Fixed row ownership, ms |
|---|---:|---:|---:|
| Olympus z4 | 36.377 | 30.436 | 32.932 |
| Norne z8 | 20.006 | 15.623 | 24.221 |

All variants require five pressure GMRES iterations and produce identical true
relative residuals (0.000852596 / 0.00147409). Raw arrays reduce cycle time by
16% / 22%, with an extra 5.8 / 4.6 KB per application. Numeric update is slower
because the schedule currently validates copied graph buffers in the prototype.
Fixed row ownership allocates about 185–187 KB/application and is not retained.
Full simulation validation of raw GS is in progress. Production would own its
schedule on the smoother state instead of using the experimental weak cache.

The application already selects native system ILU at exactly one worker. The
other session's below-four-worker change remains unmerged because the restored
worker metadata makes native ILU partition at two/three workers; that crossover
needs its own evidence. No default change is needed for the one-worker case.

The paired Float64/Int32 execution screen is complete (three fresh samples per
policy, four policies per model, alternating round order). All 24 final states
are bitwise identical within their model, with zero timed compilation and
unchanged Newton/linear counts. Medians in seconds:

| Model at 15 workers | Baseline | Typed tiles | Raw GS | Both |
|---|---:|---:|---:|---:|
| Olympus z4 | 18.923 | 18.800 | 18.093 | 18.016 |
| Norne z8 | 9.840 | 9.572 | 9.587 | 9.346 |

Combined gain is 4.8% / 5.0%, rather than the larger gain suggested by isolated
property timings. Property-update medians fall 0.492→0.334 / 0.503→0.339 s;
preconditioner application falls 10.733→9.784 / 4.364→4.113 s. No original-20%
result or OPM parity claim follows from this short diagnostic.

A small source change now selects natural CPU ILU/DILU traversal at one runtime
worker, retaining dependency levels for multiple workers. The same row loops
are reused through a function barrier for factorization, triangular application,
and fused correction. The initial one-worker correctness pass completed 68 new
checks plus the existing KA/symbolic-reuse tests. The first 15-worker test driver
incorrectly forced its natural reference through a parallel batch; the reference
now uses the existing Polyester thread-disable context. It is being rerun before
any source commit or performance conclusion.

The natural serial ILU change passes the new tests and existing KA/reuse suite:
1,155 assertions at one worker and 1,147 at 15 (36 summaries each). Three paired
fixed-matrix rounds on the larger Float64/Int32 fixtures give:

| One-worker operation | Dependency levels, ms | Natural rows, ms | Gain |
|---|---:|---:|---:|
| Olympus z4 ILU application | 64.596 | 21.011 | 3.07x |
| Olympus z4 factorization | 95.782 | 35.944 | 2.66x |
| Norne z8 ILU application | 51.107 | 29.588 | 1.73x |
| Norne z8 factorization | 95.298 | 69.860 | 1.36x |

All factor arrays, inverse diagonals and solution vectors agree exactly in
these probes; all timed kernels allocate zero bytes and report zero compilation.
At 15 workers both schedule selections use the same level arrays and helper:
application medians are 10.693/10.544 ms for Olympus and 10.959/11.579 for Norne;
factorization is 23.899/23.517 and 30.751/30.798 ms. These paired selections do
not measure different algorithms at 15 workers. Source was the uncommitted
natural-traversal patch on c1379eaf, retained with SHA256 in the raw results.
The application already uses faster native serial ILU by default, so this
improves explicit KA serial use without claiming a fastest-serial full-run gain.

## OPM Float64 reference preparation

Local software is Flow 2026.04 (`libopm-simulators-bin 2026.04-1~noble`) and
Open MPI 4.1.6. `prepare_opm.py` reuses OPM's parser for the two targeted decks
and mirrors the existing Julia refinement. Both exported decks pass OPM parsing
and the independent GeoEnergyIO cross-check, covering all property/PVT arrays,
faults, MINPV, completions, well/group controls and dates. The check includes
314 Olympus and 211 Norne completion records in the selected schedule prefix.
Norne's fifty reports are the ceiling of 20% of its 247 source reports.

Input-deck integer/property arrays round-trip exactly through OPM. Interpolated
ZCORN decimal round-trip differences are at most 2.221e-16 / 2.034e-16 relative.
The Julia helper inserts default MINPV before dividing it among refined layers;
the exporter now makes this value explicit for Olympus too. It also expands
saturation/PVT/depth table rows and preserves PVT region terminators, rather
than relying on a flattened OPM serialization accepted only by OPM.

Flow dry runs report exactly the same active-cell counts as Jutul:
771,000 Olympus / 355,444 Norne. Norne removes 3,972 of the 359,416 input-active
cells through MINPV, matching the Julia mesh. Input pore-volume totals reported
by Flow are 219,856,814 / 673,248,728 m3. Unsupported summary requests produce
warnings but do not prevent setup. No simulation comparison has been claimed.

Flow normally permits Float32 linear solves up to 20 days and uses CNV 0.01.
The reference disables single precision with `--max-single-precision-days=0`
and sets CNV 0.001, MB 1e-7 and linear reduction 0.005, including relaxed values,
to match the Julia harness. Its existing solver/timestep policies are retained;
Flow reports that it ignores SCHEDULE TUNING by default. Numerical iteration
and ministep counts will be reported rather than silently treating controllers
as identical. Binary output is disabled in the timed reference. Full wall time
and Flow's solver time must be distinguished from Julia's warmed simulation timer.

The first 15-rank reference is running, pinned to the same physical cores as the
Julia tests, one OpenMP thread per rank. Setup used about 10 GiB at the first
check; the fifteen rank CPU IDs were 0,2,...,28. Each command and deck SHA256 is
retained with Flow's version, stdout, process timing and PRT/INFOSTEP diagnostics.

The first actual 15-rank Olympus Flow diagnostic completed with process exit
zero in 681.60 s wall time (11:21.60), including input/setup. Its stdout option
suppressed the final solver-work summary, so this is not yet a solver-time or
iteration-normalized comparison with Jutul's 857.859 s simulation timer. The
controller failed *after* Flow returned because the runner file had been edited
while its shell was waiting; reading from the old file offset entered a changed
command tail. The numerical timing/command are retained. The runner now parses
its whole body as a compound block and exits from that block. Statistics are
enabled in redirected output for subsequent runs. The interrupted queue did not
start Norne or either one-rank run.

A 30-line library trial now implements the validated typed property tiling,
using the existing minbatch and worker interface. It is restricted to large
batch-threaded KA CPU updates with common entity counts; original paths handle
serial, small and mixed-entity updates. Focused dependency/AD/matrix-layout and
one/15-worker context tests, followed by both three-sample Float64/Int32 short
simulations, are running before a production commit. Startup source commits and
source-diff SHA256 are now captured by the harness. Final active-cell maps are
saved outside the timer to support correctly aligned OPM field comparisons.

A separate waiting controller resumes OPM only if both production-trial complete
state comparisons are bitwise equal. It runs Norne MPI15, separate double-output
MPI15 state validation for both cases, then the one-rank references. These
remain sequential; no numerical benchmark overlaps another numerical run.

The property-tiling trial passed the focused and existing context/reservoir
checks at both worker counts: 603 / 689 assertions across 32 summaries each.
Its new test contributes 30 / 90 assertions, covering dependent scalar/matrix
properties, Float64 forward derivatives, small/mixed-entity fallbacks, actual
tile dispatch, and context worker limits. The existing SPE1 KA adjoint checks
also pass. Full refined state/performance validation is now warming; no property
source commit has been made yet. Tiling stays inside each worker's existing
batch, so it does not launch a thread batch for each cache tile.

Committed property tiling in Jutul `ddeb7da9` after the full validation. It adds
30 library lines and no public control knobs, using one existing worker batch
and multiple local cache tiles. The two three-repeat Float64/Int32 short runs
finish at medians 18.359 s Olympus / 9.231 s Norne, with zero timed compilation,
25/157 and 17/82 Newton/linear iterations, no failed work, and bitwise-identical
complete pressure/saturation arrays. Secondary-update medians are 0.328 / 0.337 s,
about one third below the earlier paired untiled 0.492 / 0.503 s. Treat whole-run
comparisons across batches cautiously; the earlier paired tiling-only gain was
0.65% / 2.73%. These runs demonstrate the implementation and state preservation,
not final fivefold scaling or OPM parity.

The superseded property override is removed from the performance harness.
`compare_execution.jl` now compares the current library baseline with raw GS;
the four-policy historical measurements remain tied to `d313001c` and their
recorded Jutul source commit. Active-cell order and phase-name metadata were
saved successfully for both cases. The OPM continuation gate passed and its
Norne 15-rank reference is running with solver diagnostics redirected to files.

## First complete Norne MPI15 Float64 reference

Flow 2026.04 completed the matched refined Norne first-50-report schedule:
846.95 s simulation, 11.72 s setup (4.14 s deck input), 859.43 s process wall.
It uses 15 pinned ranks, one thread each, Float64 linear solves and the recorded
common CNV/MB/linear tolerances. Solver statistics are now visible in redirected
stdout. There are 267 timesteps, 2,620 linearizations (630 wasted), 2,383 Newton
iterations (630 wasted), and 2,097 linear iterations (359 wasted).

| Norne original 20% phase | Flow MPI15, s | Earlier Jutul KA15, s |
|---|---:|---:|
| Simulation | 846.95 | 1171.912 |
| Assembly | 253.78 | 179.12 |
| Properties/update | 69.58 | secondary 83.86 |
| Linear solve, including setup | 500.96 | setup+preconditioner+Krylov about 807.25 |
| Linear setup (subset) | 271.59 | 207.64 |

Phase definitions and Krylov work units differ; keep these counts alongside
timing. Existing nonlinear/timestep policies remain in place. This Jutul
reference predates the new cache tiles and uses Int64 indices. It cannot be
substituted for a new matched production run or a one/15-worker scaling pair.
No fivefold or parity claim is made.

Flow's actual logged pressure configuration has `add_wells=true` and
`use_well_weights=false`. The earlier Jutul reference uses `:cpr`, while the
existing `:cprw` option includes well pressure coupling. A single targeted
Float64 comparison now checks this existing interface with GMRES, KA ILU0,
pressure GS and Int32 indices fixed. There is no new solver implementation or
parameter sweep. The short simulations are bounded and will retain full-state
comparisons before any longer validation. The rest of the MPI controller is
held at the completed-run boundary, so measurements remain sequential.

The targeted existing `:cprw` Float64 option is rejected for these fixtures.
Olympus median is 18.395 s versus the current library's 18.359 s, with the same
25 Newton iterations but 159 versus 157 linear iterations and about 400 MB more
allocation. Its final maximum differences are 192.55 Pa / 7.06e-6 saturation.
Norne regresses to median 41.147 s versus 9.231 s, with 98 / 311 Newton/linear
iterations (63 / 190 wasted), compared with 17 / 82 and no waste. Final maximum
differences reach 24,741 Pa / 0.004984 saturation after the different solve
history. All timed samples have zero compilation; none justifies promoting the
variant. The application default and nonlinear/timestep settings are unchanged.
No further configuration sweep is planned from this result.

Resumed the separate Olympus MPI15 output-validation run with double-precision
binary output, retaining its existing Flow solver and strict reference tolerances.
The controller will hold again after this run. The existing OPM binary reader
and saved Jutul active-cell ordering will support the complete-field check.

## Olympus MPI15 fields and Float64 follow-up

The separate Olympus output-validation run completed in 667.26 s simulation,
19.44 s setup and 688.27 s process wall, with 148 timesteps. It reports 1,640
linearizations (777 wasted), 1,529 Newton iterations (777 wasted) and 1,414
linear iterations (534 wasted). Assembly is 142.40 s, linear solve 422.71 s
(including 196.79 s setup), properties/update 79.65 s and output 3.31 s. These
are Flow's phase definitions and work units. This is a second reference run,
with output enabled; it does not replace the output-disabled wall sample.

An official `opm.io` restart extraction checks metric units, report 4 and day
1461, and matches all 771,000 global active-cell identities before permuting
fields into the saved Jutul mesh order. Against the earlier original-20% Jutul
state, relative pressure norm difference is 1.52323e-4 and maximum difference
246,892 Pa. Saturation relative norm difference is 5.11604e-4, maximum 0.0587591.
This establishes a complete aligned comparison, not numerical equivalence.
The restart pressure and saturation arrays are Float32 despite the requested
`--ecl-output-double-precision=true`; this is recorded separately from Flow's
Float64 linear algebra (`--max-single-precision-days=0`). No timestep or nonlinear
control experiments follow from these field differences.

The Norne MPI15 output-validation run is next, with the parent controller held
at its following boundary to keep all numerical jobs sequential. A focused
Float64 numeric-ILU trial uses the already stored diagonal positions to skip
irrelevant entries during factorization. The previously rejected triangle-bound
*solve* optimization is excluded; no repeat of that configuration is needed.
The factor-only patch awaits paired measurements and correctness checks.
Remote `sol-perf` was checked again: Jutul c1379eaf / Darcy 2d5ee982 unchanged.

A second fixed-array probe targets black-oil increment diagnostics: the current
CPU implementation reads the same increment/state arrays six times for three
phase sums and maxima. A fused Float64 host trial is confined to `results/`,
with empty/strided/phase/NaN/Inf checks and matching exact maxima; it is not a
nonlinear or timestep-control change. Queue order is OPM fields, ILU-factor
replay, then increment replay. The original OPM parent remains held until these
short investigations finish, so no one-rank reference overlaps Julia work.

The factor-only diagonal-bound trial is rejected and its library patch reverted.
Three paired fixed-matrix rounds preserve factors/diagonals/solutions bitwise,
with zero timed compilation/allocation. Serial medians improve 37.208→33.822 ms
Olympus and 71.456→66.027 ms Norne, but the active 15-worker path is 22.933→23.114
ms / 30.718→30.122 ms. This does not justify another library change for the
current performance target. Raw results and the rejected patch remain in results.

The increment diagnostic probe passes 624 checks at each worker count. On
355,444 strided Float64 increments with dual-valued black-oil state, the original
six reductions allocate 17,062,000 bytes per call. Medians are 4.811→1.195 ms at
one worker and 5.899→0.085900 ms at 15. The trial reports 64 bytes per call and
zero timed compilation. A guarded 30-line library implementation now awaits
focused/layout/backend-fallback tests and complete one/15-worker state runs.
Maxima are compared exactly; sums are compared with floating-point tolerance.
No timestep or nonlinear controls are changed.

Norne's output-validation run completes in 854.48 s simulation, 12.33 s setup,
867.77 s process wall; output is 6.14 s. Work counts match the output-disabled
run exactly: 2,620 linearizations, 2,383 Newton, 2,097 linear iterations, 267
steps. All 355,444 active-cell identities are present, but the latest restart
is report 44 / day 541, while the simulation finishes report 50 / day 572.
The field comparator correctly rejects that endpoint. A final-restart output
request is needed before making a Norne cross-engine field comparison.

The first guarded increment implementation passes 8,913 / 8,923 assertions
(19 summaries at one / 15 workers), including 8,643 / 8,645 focused checks and
the existing KA adjoint simulation. Its short full-run medians are 20.079 / 9.424 s,
17 Newton / 82 linear iterations, no wasted work or timed compilation. The
15-worker complete state is bitwise equal; the serial state differs from the
parallel reference only by 4.01e-6 Pa / 3.69e-12 saturation. However, allocations
remain about 2.765 / 2.705 GB: this reveals that the real simulator misses its
`StridedVector` guard. These timings do not validate a production speedup.

The actual cell-major view is a one-row `SubArray` selected with a range, whose
`vec` is a noncontiguous `ReshapedArray`. A new regression reproduces the missed
path: 8,645 checks pass and its allocation assertion fails (960,656 bytes at
20,000 cells). The guard now also recognizes a reshaped vector's strided host
parent; the reduction uses checked dimensions and explicit linear indices.
Focused/convergence one/15-worker tests and fresh Norne simulations are being
repeated under the distinct `production_increment_views` label. Both OPM
controllers are held while this correction is validated; no numerical overlap.

The harness's existing solver-option forwarding is independently validated:
three Float64/Int32 Olympus z4 samples with the previously measured four-sweep
SPAI configuration give median 15.876 s, 23 Newton / 135 linear iterations,
zero wasted work and zero timed compilation. Complete pressure/saturation fields
are bitwise equal to that configuration's earlier saved state. `main` now accepts
`solver_options`, passes them to warm/measured/profile simulator construction,
and records their representation and actual smoother metadata. No new library
or environment option is added.

## Retained Float64 increment reduction

The final implementation uses the existing CPU linearization interface to
return a strided row view for a one-row Float64 matrix. Passing the generic
reshaped wrapper directly to Polyester had failed its 15-worker regression
in StrideArraysCore; the row view avoids that unsupported representation.
The library change is 34 lines, with no state fields, public option or new
solver. Large Float64 host increments with the full configured worker count
use one six-accumulator batch; small, other-precision/backend and limited-worker
inputs retain the existing reductions. Dimensions are checked before batching.

The actual cell-major regression now compares independently against `vec`
values, and passes its allocation assertion. Focused checks pass 8,646 / 8,648
at one / 15 workers, plus 46 / 54 convergence checks. A final current-source
15-worker KA run also passes its reservoir, precision/backend, hybrid and
26 adjoint assertions. Complete short Norne states are bitwise equal to the
immediately preceding same-index baseline at both one and 15 workers. The older
Int64 serial reference differs only by 9.38e-6 Pa / 7.69e-12 saturation.

Fresh three-sample runs retain 17 Newton / 82 linear iterations, no waste,
and zero timed compilation. Allocation falls 2.765→2.475 GB at one worker and
2.705→2.415 GB at 15, about 290 MB per simulation. Serial median is 20.063 s,
essentially unchanged. Its state-update median falls 0.328→0.257 s; threaded
updates fall 0.198→0.118 s. Across-batch threaded whole-run medians are noisy
(9.424 / 9.641 s before/after), so they do not establish an overall gain.

A same-process paired comparison alternates the original six-pass/vec path and
the final fused path for three rounds, with fresh simulators and BLAS=1. All
six complete states are bitwise identical, work counts agree, and compilation
is zero. Whole-run medians are 9.40435→9.34935 s, a modest 0.6% gain. Update
medians are 0.210633→0.118894 s, about 44% lower; allocations fall by about
290 MB. This justifies retaining the small allocation/update improvement,
without attributing a large application speedup to the isolated reduction.
The paired driver is a process-local adaptation of `compare_execution.jl`;
it is kept with raw results and does not alter library defaults.
