#!/bin/bash
#SBATCH --job-name=uncertain-build
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=01:15:00

# Run separately for each exported source snapshot; use fresh SLURM logs.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
module load Python/3.12.3-GCCcore-13.3.0 PyYAML/6.0.2-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-uncertain-products-20260927/build_snapshot.py build "$@"
