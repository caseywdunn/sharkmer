#!/bin/bash
#SBATCH --job-name=v4-assembly-access
#SBATCH --partition=devel
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=tmp/v4-assembly-access-20261110/logs/checks-%j.out
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
git ls-files -z src Cargo.toml Cargo.lock scripts/release_check.py | xargs -0 sha256sum
sha256sum src/kmer/access.rs src/kmer/primers.rs
scripts/release-check --offline --only fmt clippy rust python references panels fixtures
date --iso-8601=seconds
git ls-files -z src Cargo.toml Cargo.lock scripts/release_check.py | xargs -0 sha256sum
sha256sum src/kmer/access.rs src/kmer/primers.rs "$CARGO_TARGET_DIR/release/sharkmer"
