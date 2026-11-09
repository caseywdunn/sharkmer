#!/bin/bash
#SBATCH --job-name=uncertain-compare
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=05:30:00

# Arguments are forwarded to run_allocated_comparison.py. This job requires
# completed compiled checks, immutable build receipts and final source review.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
module load Python/3.12.3-GCCcore-13.3.0 PyYAML/6.0.2-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-uncertain-products-20260927/run_allocated_comparison.py "$@"
