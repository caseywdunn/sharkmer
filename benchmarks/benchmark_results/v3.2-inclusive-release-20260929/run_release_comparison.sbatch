#!/bin/bash
#SBATCH --job-name=inclusive-release
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=08:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module load Python/3.12.3-GCCcore-13.3.0 PyYAML/6.0.2-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-inclusive-release-20260929/run_release_comparison.py "$@"
