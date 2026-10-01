#!/bin/bash
#SBATCH --job-name=sequence-retention-replay
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=08:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module load Python/3.12.3-GCCcore-13.3.0 PyYAML/6.0.2-GCCcore-13.3.0
exec python3 benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/run_comparison.py \
    --old-binary tmp/inclusive-release-20260929/build-v1/candidate/sharkmer \
    --new-binary tmp/sequence-retention-20261001/checks-v2/sharkmer "$@"
