#!/bin/bash
#SBATCH --job-name=sharkmer-v4-baseline-verify
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=tmp/v4-performance-baseline-20261109/logs/verify-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
python benchmarks/benchmark_results/v4-performance-baseline-20261109/verify_baseline.py
