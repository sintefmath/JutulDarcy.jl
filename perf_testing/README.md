# KA CPU thread scalability

Working log for the `sol-perf` branches of Jutul and JutulDarcy. All measurements
use the original checkouts (`~/juliadev/Jutul` and `~/juliadev/JutulDarcy`), not
Jutul2/JutulDarcy2. Existing example edits are excluded from these commits.

See the numbered Markdown files for experiments, results, and decisions.
Benchmark scripts and a small standalone Julia environment will live here.
Generated case caches and raw profiles should not be committed.

The main fixed-SPAI0 results are 3.22x scaling for Olympus z4 (771,000 cells)
and 3.01x for Norne z8 (355,444 cells), from one to ten Julia workers. The clean
whole-simulation 4x target was not reached. Ten workers is better than fourteen
on this M4 Pro. See 006–009 for the complete sweeps and solver-policy effects.

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
All runs use `mode=:ka`, GMRES, and Float64. BLAS defaults to one thread;
`PERF_BLAS_THREADS` supports a separate dense-vector threading experiment.

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

## Choosing the pressure smoother

The controlled sweep fixes pressure SPAI0 at every worker count. The application
instead defaults to GS at one CPU worker and SPAI0 at multiple workers. These
are different algorithms, so default-to-default timing can differ substantially
from the controlled scaling table.

For the short refined Norne fixture, GS avoids recovered cuts and is faster at
ten workers: 5.215 s versus 13.024 s with SPAI0. For Olympus z4, SPAI0 is faster:
7.143 s versus 10.175 s with serial GS sweeps. Reproduce the favored Norne setup
without changing the library's general pressure-smoother default:

```sh
PERF_PRESSURE_SMOOTHER=gauss_seidel PERF_THREADS="1 2 4 6 8 10" bash perf_testing/run_scaling.sh norne 8 norne_gs
```

Keep BLAS=1. Ten BLAS workers regressed a full Olympus sample from roughly 7 to
68 s despite faster isolated vector operations. The system ILU choice remains
automatic: serial below four workers or 100,000 cells, level-scheduled otherwise,
with explicit smoother/options taking precedence. The supplied environment and
scripts are performance fixtures; the layer refinement is scoped to these decks.
