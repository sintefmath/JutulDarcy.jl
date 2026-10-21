# 009 — Remaining bottlenecks and rejected micro-optimizations (2026-10-21)

## Pressure smoothing at ten workers

Same three-report fixtures, one BLAS thread, native parallel system ILU.
Medians of three fresh warmed simulations.

| Case | Serial GS, s | Ten-worker GS, s | Ten-worker SPAI0, s |
|---|---:|---:|---:|
| Olympus z4 | 19.464 | 10.175 | 7.143 |
| Norne z8 | 11.899 | 5.215 | 13.024 |

GS at ten workers preserves the serial GS iteration and timestep-cut counts:
Olympus 25 Newton / 157 linear; Norne 17 / 82, without recovered cuts. Norne's
absolute runtime improves 2.50x over threaded SPAI0 by avoiding its additional
iterations and recovered cuts. This is an algorithm/configuration improvement,
not a 2.50x thread-scaling gain. GS itself still performs serial sweeps.
Olympus shows why changing every threaded CPU solve to serial GS would be a
poor default: it costs about 42% more than SPAI0 here.

Investigate exact dependency levels for GS sweeps next. A correct schedule must
preserve the original row order on every structural edge, including asymmetric
patterns. Considering only lower outgoing edges is insufficient: an upper-only
edge also prevents its two rows from being evaluated concurrently. Both forward
and reverse sweeps and correction application need this guarantee. A process-
local prototype and directed-pattern comparisons are queued before any library
change.

## Full-run BLAS threading: rejected

Ten BLAS threads speed up isolated dense dot/axpy operations but make the full
Olympus z4 run dramatically slower. After warm-up, the first measured sample
takes 68.097 s, with 62.399 s in preconditioner application, versus about 7.143 s
total / 3.013 s preconditioning with BLAS=1. Iteration counts agree and measured
allocations remain around 4.5 GB. Stop the remaining repeats: the regression
is far beyond measurement noise. Mixing dense-solver/BLAS workers with frequent
Polyester operations is a likely contributor; do not change global BLAS settings
in the library. Keep benchmark BLAS=1.

## CPR row loops: rejected

Repeated operations on the same Olympus z4 matrix, ten workers:

| Operation | Existing KA CPU, ms | Host batches, ms |
|---|---:|---:|
| Pressure matrix reduction | 1.159 | 1.163 |
| True-IMPES weights | 0.230 | 0.224 |

There is no useful measured gain to justify another host implementation.

## ILU triangle bounds

Restricting CPU ILU solves to the stored lower/upper triangles is mathematically
valid because symbolic setup requires sorted columns. It improves isolated
application from 5.224 ms to 5.008 ms; restoring the original gives 5.262 ms.
However, the full Olympus z4 samples remain about 7.14 s. Do not claim a full-run
win from this microbenchmark. Retain only if the remaining full comparisons
show a useful gain.

## Largest fixed-SPAI0 fixture

Olympus z8, 1,542,000 cells: serial median 54.594 s (50.294–61.973 s), ten-worker
median 17.279 s, or 3.16x. Counts agree at 31 Newton / 251 linear, without cuts.
The wide serial spread and observed active memory compression make this less
clean evidence than z4; it does not demonstrate the 4x target.
