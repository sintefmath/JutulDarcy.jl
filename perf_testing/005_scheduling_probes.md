# 005 — Scheduling probes and rejected changes (2026-10-21)

Use repeated warmed operations inside a single process to distinguish policy
cost from fresh-model setup and process-to-process variation. Each timing below
is a median of three measurements after warm-up. Runtime method replacements
were confined to the probe process; repository source was not changed.

## AMG row and setup scheduling

Olympus z2, 385,500 cells, 10 workers. Rebuild refreshes the hierarchy in `:memory`
mode; application is one pressure preconditioner cycle.

| Worker QoS | Numeric rows | Scratch-row setup | Rebuild, ms | Apply, ms |
|---|---|---|---:|---:|
| original | batch | batch | 117.20 | 1.860 |
| original | threads | threads | 106.74 | 2.180 |
| original | batch | threads | 145.13 | 1.860 |
| interactive | batch | batch | 114.40 | 1.860 |
| interactive | threads | threads | 107.64 | 2.228 |
| interactive | batch | threads | 145.62 | 1.859 |

Polyester numeric rows improve application over Julia tasks. Static Julia tasks
can make setup a little faster in isolation, but mixing them with Polyester
numeric rows is slower here. Keep the consistent Polyester policy. Promoting
workers from default to interactive macOS QoS does not materially change the
result, so do not add process-priority changes to the benchmark or library.

## ILU batch size

The original schedule has 274 lower-dependency levels on Olympus z2 and 209 on
Norne z2. Median level widths are 1,426.5 and 383 rows, respectively. Increase
Polyester's minimum batch from 1 to 16/32/64/128/256 and compare the same fixed
factorization and RHS. Full solution arrays match bit for bit across variants.

| Minimum rows per worker | Olympus apply, ms | Norne apply, ms |
|---:|---:|---:|
| 1 | 2.646 | 1.084 |
| 16 | 2.619 | 1.067 |
| 32 | 2.702 | 1.094 |
| 64 | 2.749 | 1.326 |
| 128 | 2.991 | 2.301 |
| 256 | 3.940 | 4.634 |

Larger batches do not help these cases and can substantially regress Norne.
Keep the existing fine ILU batching; do not add a tuning constant based on an
assumed universal launch-cost threshold.

## CPR CPU reduction and default system smoother

JutulDarcy commit `aa2fe00e` selects level-scheduled ILU0 automatically for KA
CPU models with at least 100,000 cells and multiple workers, when the system
smoother/options have not been explicitly supplied. Smaller cases and explicit
options retain their existing choice. This threshold brackets the measured
88,861-cell Norne case (small ILU difference) and 385,500-cell Olympus case
(clear benefit); the thread-count sweep will further check the policy.

The pressure RHS now uses a host Polyester reduction with the component count
specialized, sharing its row calculation with the accelerator kernel. On the
same forced-parallel-ILU refined fixtures, the 10-worker medians improve from
2.653 to 2.547 s for Olympus, and from 0.922 to 0.827 s for Norne. Newton/linear
iteration counts remain unchanged. The existing complete reservoir KA test
file passes at 10 workers, including hybrid models, compositional wells, and
the SPE1 adjoint solve.

Next: compare 1 and 10 workers on larger Olympus z4 and Norne z8 grids, then
measure intermediate worker counts and 14 workers.
