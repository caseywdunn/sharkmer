#!/bin/bash
#SBATCH --job-name=sharkmer-restore10
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=tmp/inclusive-release-20260929/restore-%j.slurm.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module load Python/3.12.3-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/restore_inputs.py "$@"
