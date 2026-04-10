#!/usr/bin/env bash
# Submit SLURM jobs for sharkmer panel validation on Bouchet.
#
# Two modes:
#
#   Default (one job per panel, each runs the panel's full sample x depth set):
#     bash scripts/validate_all_panels_slurm.sh [--no-blast]
#
#   Sweep mode (one job per panel x knob x value cell, each runs the panel's
#   samples at the highest declared depth tier with the swept knob applied):
#     bash scripts/validate_all_panels_slurm.sh --sweep
#
# The sweep matrix is hardcoded below in the SWEEPS array. Each entry is
# "knob:value1 value2 ...". The script submits one job per (panel, knob,
# value) combination, so the total job count is roughly:
#   n_panels x sum(values per knob)
# For the default Tier 1 matrix (4 knobs x 4 values) across 8 panels that
# is 128 jobs.
#
# Use tmux so job monitoring survives SSH disconnects:
#   tmux new -s validate
#   bash scripts/validate_all_panels_slurm.sh --sweep
#   # Ctrl-b d to detach; tmux attach -t validate to reconnect
#
# Prerequisites:
#   - sharkmer-bench conda environment exists
#   - sharkmer binary is built (cargo build --release)
#   - logs/ directory is created (script does this automatically)
#
# Build note (Bouchet):
#   module load GCC/13.3.0 Rust
#   cargo build --release
# GCC/13.3.0 provides a modern assembler (--gdwarf-5 support) needed
# by the ring crate. The default system assembler is too old.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PANELS_DIR="$REPO_ROOT/panels"
SCRIPT="$REPO_ROOT/scripts/validate_panel.py"
LOGS_DIR="$REPO_ROOT/logs"

# -----------------------------------------------------------------------------
# Tier 1 sweep matrix
#
# Each entry is "knob_name:value1 value2 value3 value4". The knob name maps
# to the corresponding sharkmer CLI flag: e.g. "k" -> "-k", and any other
# knob -> "--<knob>". Values are space-separated. To add or remove knobs,
# edit this array directly.
#
# When adding a new knob: make sure sharkmer actually accepts it as a CLI
# argument (`sharkmer --help-pcr` lists hidden PCR knobs) and that the
# default value is reflected in the value list so we can compare against it.
# -----------------------------------------------------------------------------
SWEEPS=(
    "k:17 19 21 23"
    "max-primer-kmers:10 20 40 80"
    "high-coverage-ratio:3.0 5.0 10.0 20.0"
    "tip-coverage-fraction:0.05 0.1 0.2 0.4"
)

# -----------------------------------------------------------------------------
# Mode parsing
# -----------------------------------------------------------------------------
SWEEP_MODE=0
PASSTHROUGH_ARGS=()
for arg in "$@"; do
    if [ "$arg" = "--sweep" ]; then
        SWEEP_MODE=1
    else
        PASSTHROUGH_ARGS+=("$arg")
    fi
done

panels=("$PANELS_DIR"/*.yaml)
if [ ${#panels[@]} -eq 0 ]; then
    echo "No panel YAML files found in $PANELS_DIR"
    exit 1
fi

mkdir -p "$LOGS_DIR"

# -----------------------------------------------------------------------------
# Helper: submit one SLURM job that runs validate_panel.py with the given
# panel and extra command-line arguments. The job script is written to a
# tempfile so the heredoc interpolates at submission time, not at job
# execution time.
# -----------------------------------------------------------------------------
submit_job() {
    local panel_path="$1"
    local job_name="$2"
    local validate_args="$3"

    local tmpscript
    tmpscript=$(mktemp "/tmp/val_${job_name}.XXXXXX.sh")
    cat > "$tmpscript" <<EOF
#!/bin/bash
module reset
module load miniconda
conda activate sharkmer-bench
python $SCRIPT $panel_path $validate_args
EOF
    sbatch \
        --job-name="$job_name" \
        --partition=day \
        --time=4:00:00 \
        --cpus-per-task=8 \
        --mem-per-cpu=5G \
        --output="$LOGS_DIR/${job_name}_%j.out" \
        "$tmpscript" >/dev/null
    echo "  Submitted: $job_name"
}

if [ "$SWEEP_MODE" -eq 0 ]; then
    # -------------------------------------------------------------------------
    # Default mode: one job per panel, run all sample x depth combos.
    # Any extra args (e.g. --no-blast) are passed straight through.
    # -------------------------------------------------------------------------
    echo "Submitting ${#panels[@]} validation jobs to Bouchet..."
    echo ""
    for panel in "${panels[@]}"; do
        name="$(basename "$panel" .yaml)"
        submit_job "$panel" "val_${name}" "${PASSTHROUGH_ARGS[*]:-}"
    done
else
    # -------------------------------------------------------------------------
    # Sweep mode: one job per (panel, knob, value) cell. Each cell runs the
    # panel's samples at the highest declared depth tier only, with the
    # swept knob added via --extra-args (or via -k for the kmer length).
    # BLAST against the panel's reference DB stays on by default — it is
    # fast against the small in-panel reference DB and gives the on-target
    # signal that the sweep summary uses to pick winners.
    # -------------------------------------------------------------------------
    n_jobs=0
    for sweep in "${SWEEPS[@]}"; do
        knob="${sweep%%:*}"
        values="${sweep#*:}"
        for value in $values; do
            label="sweep_${knob//-/_}_${value}"
            if [ "$knob" = "k" ]; then
                # `-k` is plumbed as a first-class arg in validate_panel.py,
                # not via --extra-args.
                cell_args="-k $value --label $label --max-reads-tier highest"
            else
                cell_args="--extra-args \"--$knob $value\" --label $label --max-reads-tier highest"
            fi
            for panel in "${panels[@]}"; do
                pname="$(basename "$panel" .yaml)"
                submit_job "$panel" "sw_${pname}_${knob//-/_}_${value}" "$cell_args"
                n_jobs=$((n_jobs + 1))
            done
        done
    done
    echo ""
    echo "Submitted $n_jobs sweep jobs across ${#panels[@]} panels and ${#SWEEPS[@]} knobs."
fi

echo ""
echo "Monitor with: squeue --me"
echo "Logs in: $LOGS_DIR/"
