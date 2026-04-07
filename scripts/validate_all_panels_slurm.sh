#!/usr/bin/env bash
# Submit one SLURM job per panel YAML file on Bouchet.
# Each job runs validate_panel.py for a single panel.
#
# Usage (from repo root, on Bouchet login node):
#   conda activate sharkmer-bench
#   bash scripts/validate_all_panels_slurm.sh [--no-blast]
#
# Use tmux so job monitoring survives SSH disconnects:
#   tmux new -s validate
#   bash scripts/validate_all_panels_slurm.sh [--no-blast]
#   bash scripts/validate_all_panels_slurm.sh
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

# Pass through any arguments (e.g. --no-blast) to validate_panel.py.
EXTRA_ARGS=("$@")

panels=("$PANELS_DIR"/*.yaml)

if [ ${#panels[@]} -eq 0 ]; then
    echo "No panel YAML files found in $PANELS_DIR"
    exit 1
fi

mkdir -p "$LOGS_DIR"

echo "Submitting ${#panels[@]} validation jobs to Bouchet..."
echo ""

for panel in "${panels[@]}"; do
    name="$(basename "$panel" .yaml)"
    tmpscript=$(mktemp /tmp/val_${name}.XXXXXX.sh)
    cat > "$tmpscript" <<EOF
#!/bin/bash
module reset
module load miniconda
conda activate sharkmer-bench
python $SCRIPT $panel ${EXTRA_ARGS[*]:-}
EOF
    sbatch \
        --job-name="val_${name}" \
        --partition=day \
        --time=4:00:00 \
        --cpus-per-task=8 \
        --mem-per-cpu=5G \
        --output="$LOGS_DIR/validate_${name}_%j.out" \
        "$tmpscript"
    echo "  Submitted: $name"
done

echo ""
echo "Monitor with: squeue --me"
echo "Logs in: $LOGS_DIR/"
