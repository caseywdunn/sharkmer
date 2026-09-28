#!/bin/bash
#SBATCH --job-name=uncertain-probe
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=00:20:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
module load Python/3.12.3-GCCcore-13.3.0 PyYAML/6.0.2-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-uncertain-products-20260927/validation/targeted_probe.py "$@"
