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
six ministeps, and zero wasted iterations. Final state checksums match; full
array comparison remains pending. Measured 15-worker runs report zero JIT time.

The one-worker call tree identifies sparse ILU row scanning as a large cost:
both triangular solves scan the whole row and test every column despite having
stored diagonal positions and sorted columns. Probe restricted triangular ranges.
Convergence is only about 0.16–0.22 s here, so it cannot by itself close the gap.
Linux exposes sibling CPU pairs 0/1, 2/3, ..., 30/31; compare physical-core
affinity separately if scaling remains limited.

Norne z4 has 177,722 active cells. The corresponding one/15-worker medians
are 23.268 / 8.627 s (2.697x). All runs take 41 accepted Newton iterations,
143 accepted linear iterations, seven ministeps, plus 15 wasted Newton and 44
wasted linear iterations from recovered cuts. The cost of those cuts is retained.
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

Next queued workload: original-duration 20% schedules on Olympus z4 (771,000
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
bitwise-equal complete solutions and zero allocation. Fifteen-worker check pending.

The first reduction prototype captured a reassigned accumulator in a closure,
creating allocations and an unsafe reduction. Moving volume reduction to a
separate function removes that capture. An AD-array multiworker comparison
failed because the synthetic dual used a non-Jutul tag, so `Jutul.value`
correctly retained it and the reduction changed type mid-loop. Use Jutul's
existing tagged-AD constructor in the fixture and keep this out of production
until the corrected multiworker comparison passes.

## Audit status

| Area | Evidence / next check |
|---|---|
| Fused assembly and secondary properties | Existing Polyester paths; Olympus assembly scales 4.8x. Retain state comparisons. |
| Linear algebra and Krylov | Main products already use host batches. Serial BLAS/vector operations scale weakly. Other session rejected a vector wrapper due to allocation growth. |
| ILU factorization and application | Large profile share; test diagonal-bounded row traversal using existing symbolic data. |
| Pressure AMG | Host row batches already present; setup and coarse operations remain partly serial. Separate construction from repeated numeric updates in replay probes. |
| CPR weights and pressure matrix | Remaining KA CPU launches; measure host-loop dispatch without changing GPU dispatch. |
| Convergence | Multi-array Base mapreduce materializes temporaries; focused fused-reduction probe queued. |
| State updates and pressure increments | Existing strided broadcasts; measure only if profiles justify changes. |
| Wells and Schur complement | Existing per-well blocks and host/device policies; retain complete reservoir KA regression coverage. |
| Allocation, GC, output | Retain GC/allocation counters and saved-state IO in simulation timing; fixture serialization is optional and outside timed runs. |

## Concurrent work

Checked remote `sol-perf` at the user's request. Jutul is unchanged; JutulDarcy
advanced to `bc8be059`, adding M4 larger-grid results and a rejected dense-vector
wrapper probe. Its best reported scaling is about 3.22x on Olympus z4 at ten
workers. Incorporated that commit after the baseline batch. No remote writes.
