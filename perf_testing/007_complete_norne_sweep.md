# 007 — Complete Norne sweep (2026-10-21)

Production source and solver settings are the same as experiment 006.
Norne z8 has 355,444 active cells. Each sample includes three one-day report
intervals and recovered timestep cuts; compilation and input setup are excluded.

| Julia workers | Median, s | Sample range, s | Scaling |
|---:|---:|---:|---:|
| 1 | 39.232 | 39.215–39.365 | 1.00x |
| 2 | 30.990 | 30.754–31.121 | 1.27x |
| 4 | 19.362 | 19.222–20.111 | 2.03x |
| 6 | 15.773 | 15.173–15.828 | 2.49x |
| 8 | 13.706 | 13.355–14.052 | 2.86x |
| 10 | 13.024 | 12.893–13.098 | 3.01x |
| 14 | 22.656 | 22.308–23.251 | 1.73x |

Ten workers is the best measured choice for both large fixtures. Fourteen
workers slows Norne to 22.656 s (1.73x), with most of the additional cost in
preconditioning, linear setup, and other linear-solver work. This supports
using the ten performance-core worker count on this machine.

All 21 measured Norne samples retain 68 accepted Newton iterations, 269
accepted linear iterations, nine accepted ministeps, 30 wasted Newton
iterations, and 128 wasted linear iterations. Pressure/saturation checksums
are unchanged. Complete array comparisons are queued after the remaining
isolated probes.

The remaining measurements compare serial ILU at two workers, the actual
one-worker application default, CPR row loops, and threaded BLAS. No broader
Krylov vector-storage change is in either branch.
