# 004 — CPU products and parallel ILU (2026-10-21)

Jutul commit: `d24bcbcc` (pushed to `origin/sol-perf`).

Scaled KA CSR multiplication now uses a CPU row-batch implementation for
ordinary host CSR storage. It supports scalar/block values and arbitrary
alpha/beta, avoids reading the output when beta is zero, and checks dimensions.
Accelerator dispatch remains in the existing kernel/vendor path.

The existing portable `smoother_type=:ka_ilu0` provides a level-scheduled ILU0
implementation with Polyester barriers. Compared it with `:ilu0`, which selects
the existing serial CSR ILU0 in CPU KA mode. No default has been changed yet.

## Results at 10 workers

Three fresh simulations per median; three one-day reports; fixed pressure SPAI0.
All columns retain identical Newton and linear iteration counts.

| Case | Original code, s | Polyester only, s | + CPU products, s | + parallel ILU0, s |
|---|---:|---:|---:|---:|
| norne z2 | 1.092 | 1.082 | 0.936 | 0.922 |
| olympus z2 | 3.261 | 3.474 | 3.342 | 2.653 |

The scaled-product path recovers the additional Krylov cost from experiment 003.
Parallel ILU improves refined Olympus substantially; the advantage at Norne z2
is small. The 4x threading target has not yet been demonstrated. Larger grids
and the remaining setup/row scheduling costs are the next checks.

## Validation

- Scaled scalar/block products: 24 checks each at 1 and 10 workers, including
  Float32/Float64/complex values, nonzero beta, NaN output with beta zero, and
  dimension errors.
- Existing Jutul KA context/generic-AD/grouped-multimodel suite passes at
  1 and 10 workers.
- Existing Jutul KA preconditioner suite passes at 10 workers, including the
  serial-CSR versus KA ILU0 comparison, AMG, and CPU/JLArray paths.
- Complete application-state comparisons against the original code: maximum
  pressure differences below 1.1e-6 Pa, maximum saturation differences below
  9.1e-14. Iteration counts unchanged.

A scheduling diagnostic found the launching Julia thread at macOS interactive
QoS (0x21) and the other numerical workers at default QoS (0x15). A process-local
priority experiment will check whether this affects measured scaling. This is
not a core affinity guarantee and no OS-wide setting is changed.
