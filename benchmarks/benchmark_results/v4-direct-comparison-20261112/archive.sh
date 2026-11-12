#!/bin/bash
#SBATCH --job-name=v4-paired-archive
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=tmp/v4-direct-ingestion-20261110/logs/archive-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module reset
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
python tmp/v4-direct-ingestion-20261110/archive_comparisons.py
