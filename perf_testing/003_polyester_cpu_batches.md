# 003 — Polyester CPU batches (2026-10-21)

Jutul commit: `c7ea1bd3` (pushed to `origin/sol-perf`).

The generic minimum-batch calculation now uses the loop length, and the batch
implementation respects a context's requested worker count. KA CPU closure
loops use this Polyester path; derived contexts preserve the worker limit.
CPU AMG numerical row loops and scratch-row setup loops also use Polyester.
Accelerator kernels retain their existing execution path.

Validation: existing threading checks plus a small worker-limit regression
block pass with 1 and 10 default workers (23 and 49 checks, respectively).
Both refined cases retain bitwise-identical complete pressure and saturation
arrays, and unchanged Newton/linear iteration counts.

## First measured result

Three fresh runs per case, 10 workers, same ILU0/SPAI0 solver as the baseline.
Phase columns are medians of the three phase timings.

| Case | Before, s | After, s | Assembly, s | Krylov excluding preconditioner, s |
|---|---:|---:|---:|---:|
| norne z2 | 1.092 | 1.082 | 0.211 → 0.115 | 0.070 → 0.119 |
| olympus z2 | 3.261 | 3.474 | 0.316 → 0.158 | 0.451 → 0.676 |

Assembly improves substantially, but the whole refined Olympus simulation
regresses about 7%; Norne is essentially unchanged. This is an intermediate
result, not a claimed end-to-end win. Mixing frequent Polyester calls with
remaining KA CPU task launches is a hypothesis for the increased Krylov cost.
The next controlled experiment adds a host implementation for scaled CSR
products, which otherwise launch KA CPU tasks. Separately compare the existing
level-scheduled KA ILU0 against the serial system smoother.
