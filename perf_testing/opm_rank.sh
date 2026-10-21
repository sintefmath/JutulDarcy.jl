#!/usr/bin/env bash
set -euo pipefail
# Same Ryzen/WSL physical-core mapping as pin_threads.jl: sibling pairs 0/1,
# 2/3, ... . Keep one OPM rank and one OpenMP thread per physical core.
perf_rank=${OMPI_COMM_WORLD_LOCAL_RANK:-0}
perf_cpu=$((2*perf_rank))
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec taskset -c "$perf_cpu" /usr/bin/flow "$@"
