#!/bin/bash
#SBATCH --job-name=inclusive-builds
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module load Python/3.12.3-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-inclusive-release-20260929/build_pair.py "$@"
