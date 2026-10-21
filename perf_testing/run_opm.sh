#!/usr/bin/env bash
set -euo pipefail
perf_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
perf_model=${1:?case: olympus or norne}
perf_ranks=${2:-15}
perf_repeats=${3:-1}
[[ $perf_ranks -ge 1 && $perf_ranks -le 16 && $perf_repeats -ge 1 ]]
case "$perf_model" in
    olympus) perf_z=4; perf_stem=OLYMPUS_Z4 ;;
    norne) perf_z=8; perf_stem=NORNE_Z8 ;;
    *) exit 2 ;;
esac
perf_case_dir="$perf_dir/results/opm/${perf_model}_z${perf_z}"
perf_deck="$perf_case_dir/$perf_stem.DATA"
[[ -f $perf_deck && "$perf_dir/results/opm_deck_checks.toml" -nt $perf_deck ]]
for ((perf_repeat=1; perf_repeat<=perf_repeats; perf_repeat++)); do
    perf_output="$perf_case_dir/mpi${perf_ranks}_r${perf_repeat}"
    # Preserve previous results: each invocation uses fresh output directories.
    [[ ! -e $perf_output ]]
    mkdir -p "$perf_output"
    /usr/bin/flow --version > "$perf_output/version.txt"
    sha256sum "$perf_deck" > "$perf_output/input.sha256"
    perf_args=("$perf_deck" --threads-per-process=1 --max-single-precision-days=0
        --tolerance-cnv=0.001 --tolerance-cnv-relaxed=0.001
        --tolerance-mb=1e-7 --tolerance-mb-relaxed=1e-7
        --linear-solver-reduction=0.005 --relaxed-linear-solver-reduction=0.005
        --enable-ecl-output=false --enable-esmry=false --enable-terminal-output=false
        --output-extra-convergence-info=steps --output-dir="$perf_output")
    printf '%q ' mpirun -np "$perf_ranks" --bind-to none bash "$perf_dir/opm_rank.sh" "${perf_args[@]}" > "$perf_output/command.txt"
    printf '\n' >> "$perf_output/command.txt"
    /usr/bin/time -v -o "$perf_output/time.txt" \
        mpirun -np "$perf_ranks" --bind-to none bash "$perf_dir/opm_rank.sh" "${perf_args[@]}" \
        > "$perf_output/run.log" 2>&1
done
