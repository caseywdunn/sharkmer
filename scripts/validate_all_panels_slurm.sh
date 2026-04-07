#!/usr/bin/env bash
# Submit one SLURM job per panel YAML file on Bouchet.
# Each job runs validate_panel.py for a single panel.
#
# Usage (from repo root, on Bouchet login node):
#   scripts/validate_all_panels_slurm.sh [--no-blast]
#
# Prerequisites:
#   - sharkmer-bench conda environment exists
#   - sharkmer binary is built (cargo build --release)
#   - logs/ directory is created (script does this automatically)

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
    sbatch \
        --job-name="val_${name}" \
        --partition=day \
        --time=4:00:00 \
        --cpus-per-task=8 \
        --mem-per-cpu=5G \
        --output="$LOGS_DIR/validate_${name}_%j.out" \
        --wrap="
module purge
module load miniconda
conda activate sharkmer-bench
python $SCRIPT $panel ${EXTRA_ARGS[*]:-}
"
    echo "  Submitted: $name"
done

echo ""
echo "Monitor with: squeue --me"
echo "Logs in: $LOGS_DIR/"
