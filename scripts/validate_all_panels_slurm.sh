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
# Can also be submitted as a SLURM job itself (it will then submit child jobs):
#   sbatch scripts/validate_all_panels_slurm.sh [--no-blast]
# or from within the scripts/ directory:
#   sbatch validate_all_panels_slurm.sh [--no-blast]
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

# Resolve repo root robustly for both `bash scripts/…` and `sbatch …` invocations.
# When run via sbatch, SLURM copies the script to a temp path in /var/spool/slurmd/
# so dirname "$0" is wrong.  SLURM_SUBMIT_DIR is the directory where sbatch was
# called, which is either the repo root or the scripts/ subdirectory.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    if [ -f "${SLURM_SUBMIT_DIR}/Cargo.toml" ]; then
        REPO_ROOT="${SLURM_SUBMIT_DIR}"
    else
        REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    fi
else
    REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi
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
#
# Supported flags:
#   --sweep                 Sweep mode (see comments near the top).
#   --panels <a> <b> ...    Restrict to the named panels (basenames without
#                           .yaml, space-separated). Consumes all following
#                           tokens until the next flag or end of args. Useful
#                           for re-running only the panels whose max_reads
#                           matrices have changed, without resubmitting the
#                           full 128-cell sweep every time.
#   anything else           Forwarded to validate_panel.py in default mode.
SWEEP_MODE=0
PASSTHROUGH_ARGS=()
PANEL_FILTER=()
while [ $# -gt 0 ]; do
    case "$1" in
        --sweep)
            SWEEP_MODE=1
            shift
            ;;
        --panels)
            shift
            while [ $# -gt 0 ] && [[ "$1" != --* ]]; do
                PANEL_FILTER+=("$1")
                shift
            done
            ;;
        *)
            PASSTHROUGH_ARGS+=("$1")
            shift
            ;;
    esac
done

if [ ${#PANEL_FILTER[@]} -gt 0 ]; then
    panels=()
    for name in "${PANEL_FILTER[@]}"; do
        candidate="$PANELS_DIR/${name}.yaml"
        if [ ! -f "$candidate" ]; then
            echo "--panels: no such file $candidate"
            exit 1
        fi
        panels+=("$candidate")
    done
else
    panels=("$PANELS_DIR"/*.yaml)
fi
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
    local sbatch_mem="${4:-40G}"
    local sbatch_time="${5:-4:00:00}"

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
        --time="$sbatch_time" \
        --cpus-per-task=8 \
        --mem="$sbatch_mem" \
        --output="$LOGS_DIR/${job_name}_%j.out" \
        "$tmpscript" >/dev/null
    echo "  Submitted: $job_name ($sbatch_mem, $sbatch_time)"
}

# Per-panel resource overrides. Panels declaring max_reads up to 16M (see
# panels/cnidaria.yaml, panels/insecta.yaml) need ~60-80 GB for kmer counting
# at the largest depth, and the 5-depth sweep cells take appreciably longer
# than the single-depth cells. Panels not listed here fall back to the
# submit_job defaults (40G / 4:00:00).
panel_mem() {
    case "$1" in
        cnidaria|insecta) echo "96G" ;;
        *) echo "40G" ;;
    esac
}
panel_time() {
    case "$1" in
        cnidaria|insecta) echo "8:00:00" ;;
        *) echo "4:00:00" ;;
    esac
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
        submit_job "$panel" "val_${name}" "${PASSTHROUGH_ARGS[*]:-}" \
            "$(panel_mem "$name")" "$(panel_time "$name")"
    done
else
    # -------------------------------------------------------------------------
    # Sweep mode: one job per (panel, knob, value) cell. Each cell runs every
    # max_reads tier declared for the panel (so panels that want a
    # depth-stability view at knob x value = specific cell get it by
    # declaring multiple max_reads values in their YAML), with the swept
    # knob added via --extra-args (or via -k for the kmer length). BLAST
    # against the panel's reference DB stays on by default — it is fast
    # against the small in-panel reference DB and gives the on-target signal
    # that the sweep summary uses to pick winners.
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
                cell_args="-k $value --label $label --max-reads-tier all"
            else
                cell_args="--extra-args \"--$knob $value\" --label $label --max-reads-tier all"
            fi
            for panel in "${panels[@]}"; do
                pname="$(basename "$panel" .yaml)"
                submit_job "$panel" "sw_${pname}_${knob//-/_}_${value}" \
                    "$cell_args" \
                    "$(panel_mem "$pname")" "$(panel_time "$pname")"
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
