#!/bin/bash
#SBATCH --job-name=sharkmer-v4-contracts
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=tmp/v4-contracts-20261109/logs/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
module reset
source tmp/uncertain-products-20260927/toolchain/activate.sh
export CARGO_TARGET_DIR="$PWD/tmp/sequence-retention-20261001/target"
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=1 LC_ALL=C
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
export PYTHONPATH="$PWD/tmp/prerelease-checks-20261010/python-deps${PYTHONPATH:+:$PYTHONPATH}"
date --iso-8601=seconds
git rev-parse HEAD
sha256sum Cargo.toml Cargo.lock scripts/release_check.py src/kmer/identity.rs src/kmer/encoding.rs src/kmer/mod.rs
rustc --version
python --version
scripts/release-check --offline --only fmt clippy rust python references panels fixtures
date --iso-8601=seconds
sha256sum Cargo.toml Cargo.lock scripts/release_check.py src/kmer/identity.rs src/kmer/encoding.rs src/kmer/mod.rs "$CARGO_TARGET_DIR/release/sharkmer"
