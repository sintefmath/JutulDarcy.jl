# 008 — Solver policy matters as much as scheduling (2026-10-21)

These measurements use the same production source as 006/007. Three fresh
simulations per median, three one-day reports. Explicit solver overrides keep
the experiment independent of automatic policy.

## Two workers: keep serial system ILU

Pressure SPAI0 is fixed in this comparison.

| Case | Level-scheduled ILU0, s | Serial CSR ILU0, s |
|---|---:|---:|
| Olympus z4 | 17.784 | 14.984 |
| Norne z8 | 30.990 | 26.278 |

Newton, linear-iteration, and recovered-cut counts agree. The automatic
system-smoother policy should keep serial ILU at two workers. The final
policy adjustment will be recorded with the remaining probe results.

## Actual application defaults at one worker

The original application chooses Gauss–Seidel pressure smoothing at one CPU
worker and SPAI0 at multiple workers. The fixed-SPAI0 sweep deliberately
removes this algorithm change; it does not describe default-to-default scaling.

| Case | Serial default GS, s | Serial fixed SPAI0, s | Ten-worker default SPAI0, s |
|---|---:|---:|---:|
| Olympus z4 | 19.464 | 23.022 | 7.143 |
| Norne z8 | 11.899 | 39.232 | 13.024 |

Olympus GS takes 25 Newton and 157 linear iterations, versus 26/223 with SPAI0.
Norne GS takes 17 Newton and 82 linear iterations, three accepted ministeps,
and **no recovered cuts**. SPAI0 takes 68/269 plus 30 wasted Newton and 128
wasted linear iterations, with nine accepted ministeps. This makes Norne's
actual default ten-worker run slower than its actual serial default, despite
3.01x scaling when SPAI0 is fixed. Do not hide this behind the controlled table.

Next compare GS at ten workers on both fixtures. Changing the default based
only on available worker count can trade away enough preconditioner strength
to overwhelm the gain from parallel loops. Any final pressure-smoother policy
must account for both observed runtime and convergence behavior.

## Largest fixture: memory caveat

Olympus z8 has 1,542,000 active cells and takes 31 Newton / 251 linear iterations
without recovered cuts. The serial fixed-SPAI0 median is 54.594 s. Sample spread
is substantially larger than on z4; the process uses roughly 7 GB RSS and the
machine is actively compressing memory, with substantial existing swap use.
Treat this as a memory-pressure experiment, not clean evidence for a CPU-only
4x claim. The smaller z4 fixture remains the main Olympus comparison.
