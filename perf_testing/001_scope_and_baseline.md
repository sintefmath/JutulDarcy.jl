# 001 — Scope and baseline (2026-10-21)

## Starting point

- Jutul: `878698fb` on `sol-perf`.
- JutulDarcy: `a218be86` on `sol-perf`.
- Apple M4 Pro, 24 GB RAM, 10 performance cores + 4 efficiency cores.
- Julia 1.12.7; KernelAbstractions 0.9.45; Polyester 0.7.19.
- Target: at least 4x warmed simulation speedup on a sufficiently large case,
  measured at 10 default-pool Julia threads. macOS chooses core placement;
  a Julia worker count does not guarantee pinning to performance cores.

## Method

Use Norne without hysteresis and Olympus realization 1 from GeoEnergyIO artifacts.
Keep only a few report steps. Measure original grids first, then spatially refined
grids. Fix the pressure smoother across thread counts (SPAI0) to compare the same
algorithm. Also inspect the actual default solver choice separately. Keep BLAS
single-threaded. Exclude compilation and input/mesh setup; retain Newton and
linear iteration counts, convergence, allocations, and phase timings. Run CPU
benchmarks sequentially, in separate Julia processes for each worker count.

Use profiles and focused repeated-operation timings to select small changes.
Rerun the same cases after each change, checking iteration counts and final
pressure/saturation values. Add only focused correctness coverage where needed.
Commit and push incrementally with explicit `HEAD:refs/heads/sol-perf` refspecs.

## Initial code audit (hypotheses, not measurements)

- `threaded_loop_minbatch` computes `N_threads ÷ minbatch` rather than dividing
  the loop length by the batch size. With the usual batch size 1000 this chooses
  one batch, even for very large arrays.
- KA CPU closure loops use KernelAbstractions' task launches above 1000 entries.
- KA preconditioner row loops use dynamic Julia task scheduling.
- AMG setup uses static Julia task scheduling and per-thread scratch storage.
- Five-argument KA CSR multiplication always launches a KA kernel, including
  small well blocks; the three-argument CPU path already has a host specialization.

No production changes have been made yet.
