#!/bin/bash
#SBATCH --job-name=v4-direct-checks
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=tmp/v4-direct-ingestion-20261110/logs/checks-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module reset
source tmp/uncertain-products-20260927/toolchain/activate.sh
export CARGO_TARGET_DIR="$PWD/tmp/sequence-retention-20261001/target"
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=1 LC_ALL=C
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
export PYTHONPATH="$PWD/tmp/prerelease-checks-20261010/python-deps${PYTHONPATH:+:$PYTHONPATH}"
date --iso-8601=seconds
git rev-parse HEAD
rg --files -0 src tests -g '*.rs' | sort -z | xargs -0 sha256sum
sha256sum Cargo.toml Cargo.lock scripts/release_check.py
scripts/release-check --offline --only fmt clippy rust python references panels fixtures
date --iso-8601=seconds
rg --files -0 src tests -g '*.rs' | sort -z | xargs -0 sha256sum
sha256sum Cargo.toml Cargo.lock scripts/release_check.py "$CARGO_TARGET_DIR/release/sharkmer"
