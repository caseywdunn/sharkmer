#!/bin/bash
#SBATCH --job-name=sharkmer-final-checks
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
source tmp/uncertain-products-20260927/toolchain/activate.sh
export CARGO_TARGET_DIR="$PWD/tmp/sequence-retention-20261001/target"
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=1 LC_ALL=C
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
export PYTHONPATH="$PWD/tmp/prerelease-checks-20261010/python-deps${PYTHONPATH:+:$PYTHONPATH}"
exec /nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin/python benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/run_final_checks.py "$@"
