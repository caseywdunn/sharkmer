#!/bin/bash
#SBATCH --job-name=uncertain-checks
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=08:00:00

# Submit from the repository root, with --output pointing to a fresh SLURM log.
# Arguments are forwarded to run_compiled_checks.py; no comparison is run here.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
module load Python/3.12.3-GCCcore-13.3.0 PyYAML/6.0.2-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-uncertain-products-20260927/run_compiled_checks.py "$@"
