# KA CPU thread scalability

Working log for the `sol-perf` and `sol-perf-2` branches of Jutul and JutulDarcy.
Experiments 001–005 used the original macOS checkouts. Experiment 006 continues
on a Ryzen 7950X using sibling `Jutul.jl`, `JutulDarcy.jl`, and `GeoEnergyIO.jl`
checkouts in `D:/agent-jd-fixes`. Commit locally; do not push this continuation.

See the numbered Markdown files for experiments, results, and decisions.
The Ryzen continuation and current fivefold/15-worker target are documented in
[`007_ryzen_scaling.md`](007_ryzen_scaling.md).
Benchmark scripts and a small standalone Julia environment will live here.
Generated case caches and raw profiles should not be committed.

## Reproduction

Place the checkouts in sibling `Jutul.jl`, `JutulDarcy.jl`, and `GeoEnergyIO.jl`
directories. MultiComponentFlash is registry-resolved. From JutulDarcy:

```sh
julia +1.13 --project=perf_testing -e 'using Pkg; Pkg.instantiate()'
julia +1.13 --project=perf_testing --threads=1,0 --startup-file=no perf_testing/benchmark.jl olympus 2 3 3 baseline
julia +1.13 --project=perf_testing --threads=10,0 --startup-file=no perf_testing/benchmark.jl norne 2 3 3 baseline
```

Arguments are `case`, vertical refinement factor, report count, measured repeats,
and output label. `case` is `norne` or `olympus`; factor 1 is the original grid.
`PERF_MAX_DAYS` caps each report interval (default 1). `PERF_PROFILE=true`
adds one separately profiled run. `PERF_PRESSURE_SMOOTHER=default` measures the
thread-dependent application default instead of fixed SPAI0.
`PERF_SYSTEM_SMOOTHER=ilu0` selects CPU CSR ILU0 (partitioned across CPU workers
after the worker-metadata fix in the Ryzen continuation);
`PERF_SYSTEM_SMOOTHER=ka_ilu0` forces level-scheduled ILU0. The benchmark default
leaves the application choice in place.
`PERF_FRACTION=0.2` selects the first 20% of report intervals at their original
durations (rounding the report count upward), overriding the positional report
count. The default one-day cap is then disabled. Use this for final schedule
validation; the three-day runs are diagnostics. `PERF_JULIA_CHANNEL` chooses
the Juliaup channel in `run_scaling.sh` (default `1.13`).
Runs use `mode=:ka`, Float64, and GMRES by default; `PERF_MODE` and `PERF_LINEAR_SOLVER` select
another existing Krylov solver. BLAS defaults to one thread;
`PERF_BLAS_THREADS` supports a separate dense-vector threading experiment.
On Linux, `PERF_CPUSET` optionally passes an explicit CPU list to `taskset` in
the sweep runner. On this Ryzen under WSL, `0-28:2` uses 15 distinct physical
cores. Keep affinity experiments separately labeled and record the setting.

TOML results include wall/GC time, allocations, phase timings, iteration counts,
checksums, cell counts, and source commits. Serialized final pressure/saturation
arrays support complete state comparisons. No full model or dataset is committed.
Newton and linear iteration totals include failed attempts. The `wasted_*`
fields are subsets of those totals; subtract them to obtain successful-step work.

Run a thread-count sweep sequentially with:

```sh
bash perf_testing/run_scaling.sh olympus 4 scaling
PERF_THREADS="1 4 10 14" bash perf_testing/run_scaling.sh norne 8 scaling
```

The sweep defaults to 1, 2, 4, 8, and 15 workers, saving each process log
and TOML result under `results/`. Set `PERF_REPEATS` or `PERF_STEPS` to override
the default three measured samples or three report intervals. Warm-up is always
separate. The short scripts do not load plotting packages.

`compare_solvers.jl` screens the existing Krylov, system-smoother and pressure-
smoother choices. Select lists with `PERF_SOLVERS`, `PERF_SYSTEM_SMOOTHERS`,
and `PERF_PRESSURE_SMOOTHERS`. `PERF_LABEL` separates result sets.
`PERF_SCREEN_INFO=1` shows warm-up progress. Optional
`PERF_SCREEN_LIMIT_SECONDS` aborts an over-budget candidate through the existing
iteration hook and records its status before continuing; it defaults to infinity.
Use unbounded fresh runs for final performance validation.

`compare_amg.jl` uses the same driver to compare the existing AMG update APIs:
default, fixed sparsity, fixed operators, a 20% rebuild tolerance, and disabled
partial updates. `PERF_AMG_POLICIES` selects a subset by name. These are
experiments, with no change to application defaults.

`PERF_SINGLE_PRECISION_LINEAR=true` selects the existing mixed-precision solve
interface (Float64 simulation, Float32 linear algebra, FGMRES). Keep its results
separately labelled and compare complete states before selecting it.
`PERF_INDEX_TYPE=Int32` uses the existing narrower-index interface; default is
`Int`. Record this separately from floating-point precision.

For exploratory long runs, `PERF_WARMUP_MAX_DAYS=1` caps a separate warm-up
prefix (`PERF_WARMUP_STEPS=3` by default) while preserving the measured schedule.
The output records both schedules and compilation time; repeat any sample with
material JIT time. Default warm-up still uses the entire measured case.
`PERF_RUN_INFO` enables progress during diagnostic measurements; final timings
default to quiet output. `probe_pressure.jl` replays saved pressure matrices
through the existing smoothers to separate kernel cost from convergence.
Nonlinear controls and timestep selection stay at their factory defaults.
The harness fixes relaxation to false; no relaxation experiment is included.

## Optional experiments

`PERF_PIN_CPUS=0,2,...,28` pins each default Julia worker to one Linux CPU. Supply
exactly one CPU per worker; the script verifies the observed affinity. This is
separate from a process-wide `taskset` mask. `PERF_GS_LEVELS=true` enables the
isolated natural-order GS schedule prototype; `PERF_GS_GLOBAL_CHUNKS=true` also
tests fixed global row ownership. Both default off and are not library options.
Use `test_gs_levels.jl` before timing a modified prototype. `probe_pressure.jl`
replays a saved pressure matrix and RHS with selected `PERF_PRESSURE_POLICIES`;
its pressure-only GMRES measurements are screening surrogates, not full CPR or
simulation timings. `PERF_SAVE_STATES=true` saves comparison-driver final arrays
outside its measured interval.

`probe_terms.jl fixture output.toml` measures fixed-input CSR products, ILU,
CPR, Krylov vector operations and individual AMG cycle stages. Run separate
processes for each thread count, sequentially. `probe_full_terms.jl` adds the
well-reduced operator, CPR transfers and individual AMG numeric-update stages.
Isolated level timings reuse smaller working sets; their sum does not predict
a complete cycle's wall time. They do not change simulation timestepping.

`probe_builtin.jl model refinement output.toml` uses Jutul's existing
`Benchmarking.benchmark_secondary_variables`, `benchmark_assembly` and
`benchmark_linear_operators` routines on the initial physical state. Secondary
variables are timed individually, using the production host/kernel choice;
the production host path batches all variables within each worker. Retain this
distinction when comparing those timings with a full simulation.

`compare_hybrid_smoothing.jl` experiments with the existing SPAI0 and GS kernels
on different hierarchy levels. `test_hybrid_smoother.jl` checks it against a
manually assembled mixed hierarchy, including numeric and memory resets. The
prototype is confined to the performance scripts; it is not an application
default. `test_smoother_reuse.jl` invokes the retained Jutul symbolic-reuse test.
