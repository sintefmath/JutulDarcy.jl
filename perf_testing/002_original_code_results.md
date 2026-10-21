# 002 — Original code results (2026-10-21)

These are warmed medians of three fresh simulations, each with three one-day
reports. Pressure AMG uses SPAI0 at every worker count; the system smoother is
the existing CPU ILU0. BLAS uses one thread. Input and simulator construction
are outside the timer, initial AMG hierarchy construction is inside it.

| Case | Active cells | 1 worker, s | 10 workers, s | Speedup | Newton / linear iterations |
|---|---:|---:|---:|---:|---:|
| norne z1 | 44,431 | 0.807 | 0.518 | 1.56x | 9 / 31 |
| norne z2 | 88,861 | 2.356 | 1.092 | 2.16x | 15 / 52 |
| olympus z1 | 192,750 | 1.873 | 0.822 | 2.28x | 11 / 62 |
| olympus z2 | 385,500 | 7.801 | 3.261 | 2.39x | 22 / 146 |

For all four cases, complete final pressure and saturation arrays are
bitwise identical between one and ten workers, and iteration counts match. The baseline is below the 4x target.

Profiles identify these next experiments:

1. Replace the dynamic Julia task launches in AMG CPU row operations with
   Polyester, retaining the numerical operations and order within each row.
2. Use Polyester for KA CPU closure loops (including fused assembly), and fix
   the generic minimum-batch calculation and context worker limits.
3. Compare the existing level-scheduled KA ILU0 with the default serial CSR
   ILU0. The default system smoother remains serial even with 10 workers.
4. Check five-argument CSR products and AMG setup scratch-row loops for
   avoidable kernel/task launch overhead.

All layer-refined baselines are complete. Norne z2 has 88,861 active cells after
mesh processing (one fewer than twice the original active count).
Raw profiles, timings, and final arrays are retained locally in ignored
`perf_testing/results/`; numbered logs retain the reviewable conclusions.
