#!/bin/bash
#SBATCH --job-name=sharkmer-final-check
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=02:15:00
#SBATCH --output=tmp/final-check-20261108/logs/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
module reset
export LC_ALL=C
export PYTHONPATH="$PWD/scripts${PYTHONPATH:+:$PYTHONPATH}"
exec /nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin/python \
  benchmarks/benchmark_results/v3.2-final-check-20261108/run_check.py \
  --protocol benchmarks/benchmark_results/v3.2-final-check-20261108/protocol-v1.json "$@"
