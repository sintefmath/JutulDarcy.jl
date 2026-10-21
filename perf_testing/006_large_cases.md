# 006 — Larger grids and vector-operation limits (2026-10-21)

Production source: Jutul `d24bcbcc`, JutulDarcy `aa2fe00e` (the subsequent
`e9292fbd` commit only adds documentation and the sweep runner). All results use
fixed pressure SPAI0, default system-smoother selection, one BLAS thread,
three one-day report intervals, and medians of three fresh simulations after
compilation warm-up. Initial AMG construction remains inside the timer.

## Olympus: 771,000 cells (z4)

| Julia workers | Seconds | Speedup over 1 worker |
|---:|---:|---:|
| 1 | 23.022 | 1.00x |
| 2 | 17.784 | 1.29x |
| 4 | 10.489 | 2.19x |
| 6 | 8.995 | 2.56x |
| 8 | 7.490 | 3.07x |
| 10 | 7.143 | 3.22x |
| 14 | 12.731 | 1.81x |

Every run takes 26 Newton iterations, 223 linear iterations, and six accepted
ministeps, without recovered cuts. Using all 14 cores is substantially worse
than ten workers on this heterogeneous CPU. No affinity is forced; worker
count is not a guarantee of performance-core placement.

At ten workers, phase times are approximately 0.51 s assembly, 3.01 s
preconditioner application, 1.65 s other linear-solver work, and 1.28 s linear
setup. At one worker these are 3.17, 11.92, 2.75, and 3.39 s, respectively.
Assembly scales much better than the complete solver. Preconditioner barriers,
serial vector reductions, and hierarchy setup limit the overall result.

## Norne: 355,444 cells (z8), initial endpoints

| Julia workers | Seconds | Speedup over 1 worker |
|---:|---:|---:|
| 1 | 39.232 | 1.00x |
| 2 | 30.990 | 1.27x |
| 10 | 13.024 | 3.01x |

The intermediate-worker sweep is continuing. These runs have identical counts:
68 accepted Newton iterations, 269 accepted linear iterations, nine accepted
ministeps, plus 30 wasted Newton and 128 wasted linear iterations from recovered
timestep cuts. Keeping these costs in the measurement matters on this fixture.
The two-worker native ILU application is slower than the serial-worker result;
compare forced serial ILU before settling the automatic policy.

## Original-grid checks with current production code

| Case | Cells | 1 worker, s | 10 workers, s | Scaling |
|---|---:|---:|---:|---:|
| Norne z1 | 44,431 | 0.809 | 0.415 | 1.95x |
| Olympus z1 | 192,750 | 1.925 | 0.596 | 3.23x |

Compared with the original ten-worker code, these are about 20% and 27% lower
wall times. This comparison is a code improvement, not a thread-scaling ratio.
The 4x whole-simulation target has not been reached.

## Dense-vector probes: keep out of production

OpenBLAS on a 1.5-million-entry Float64 vector scales dot/axpy well, but its
norm does not: at 1/10 BLAS threads, dot takes about 0.230/0.026 ms, axpy
0.279/0.032 ms, and norm 1.436/1.448 ms. Fourteen BLAS threads also regress
dot/axpy. The harness now accepts `PERF_BLAS_THREADS` for an isolated full-run
comparison; it still defaults to one and the library changes no BLAS setting.

A private CPU vector-wrapper prototype distributes Krylov vector operations
through Polyester. Focused checks pass for real/complex operations, stable
norms, copying, and a GMRES solve. In the isolated full Norne z8 simulation it
reduces non-preconditioner solver time from about 1.54 to 0.95 s, but total
allocations rise from about 6.2 to 39.5 GB, increasing GC and secondary-property
time. Wall time regresses to 16.34 s. Runtime method replacement and associated
specialization effects complicate this experiment, but there is no validated
end-to-end gain. Do not commit the wrapper or change Krylov storage based on
its microbenchmarks. Overlapping Olympus prototype timings were discarded.

The retained production changes remain small; the broader vector-storage
experiment stays outside both repositories.
