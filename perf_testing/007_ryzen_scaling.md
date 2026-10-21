# 007 — Ryzen 7950X continuation (2026-10-21)

## Starting point and objective

Local `sol-perf-2` branches start at Jutul `d24bcbcc` and JutulDarcy
`e9292fbd`. In these Windows checkouts the SINTEF remote is named `upstream`;
the `origin` forks do not contain `sol-perf`. No pushes are authorized.
Earlier numbered results were measured on an Apple M4 Pro; do not combine
their absolute times with this workstation's results.

Target (updated by the user): at least 5x warmed simulation throughput at 15 Julia workers relative
to one worker on this Ryzen 9 7950X (16 cores / 32 hardware threads). Use refined
Olympus and Norne, with the first 20% of report intervals at their original
durations for final validation. Short three-day fixtures are diagnostic only.
BLAS stays single threaded, CPU benchmark processes run sequentially, and
warm-up is separate. Record convergence, iteration counts, phase timings,
allocations, source commits, and complete final pressure/saturation comparisons.
Compare fixed solver settings as well as the application defaults.

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
