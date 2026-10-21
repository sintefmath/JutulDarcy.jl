#!/usr/bin/env bash
set -euo pipefail
perf_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
perf_model=${1:?Usage: run_scaling.sh norne-or-olympus refinement [label]}
perf_refinement=${2:-2}
perf_label=${3:-scaling}
perf_output=${PERF_OUTPUT:-"$perf_dir/results"}
mkdir -p "$perf_output"
for perf_workers in ${PERF_THREADS:-1 2 4 8 15}; do
    perf_log="$perf_output/${perf_label}_${perf_model}_z${perf_refinement}_t${perf_workers}.log"
    julia +"${PERF_JULIA_CHANNEL:-1.13}" --project="$perf_dir" --threads="$perf_workers,0" --startup-file=no \
        "$perf_dir/benchmark.jl" "$perf_model" "$perf_refinement" \
        "${PERF_STEPS:-3}" "${PERF_REPEATS:-3}" "$perf_label" > "$perf_log" 2>&1
    tail -1 "$perf_log"
done
