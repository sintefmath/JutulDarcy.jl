# KA CPU thread scalability

Working log for the `sol-perf` branches of Jutul and JutulDarcy. All measurements
use the original checkouts (`~/juliadev/Jutul` and `~/juliadev/JutulDarcy`), not
Jutul2/JutulDarcy2. Existing example edits are excluded from these commits.

See the numbered Markdown files for experiments, results, and decisions.
Benchmark scripts and a small standalone Julia environment will live here.
Generated case caches and raw profiles should not be committed.

## Reproduction

Place the `sol-perf` checkouts in sibling `Jutul` and `JutulDarcy` directories.
The supplied environment also refers to sibling GeoEnergyIO and
MultiComponentFlash checkouts used for these measurements. From JutulDarcy:

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
`PERF_SYSTEM_SMOOTHER=ilu0` forces serial CSR ILU0;
`PERF_SYSTEM_SMOOTHER=ka_ilu0` forces level-scheduled ILU0. The benchmark default
leaves the application choice in place.
All runs use `mode=:ka`, GMRES, Float64, and single-threaded BLAS.

TOML results include wall/GC time, allocations, phase timings, iteration counts,
checksums, cell counts, and source commits. Serialized final pressure/saturation
arrays support complete state comparisons. No full model or dataset is committed.

Run a thread-count sweep sequentially with:

```sh
bash perf_testing/run_scaling.sh olympus 4 scaling
PERF_THREADS="1 4 10 14" bash perf_testing/run_scaling.sh norne 8 scaling
```

The sweep defaults to 1, 2, 4, 6, 8, 10, and 14 workers, saving each process log
and TOML result under `results/`. Set `PERF_REPEATS` or `PERF_STEPS` to override
the default three measured samples or three report intervals. Warm-up is always
separate. The short scripts do not load plotting packages.
