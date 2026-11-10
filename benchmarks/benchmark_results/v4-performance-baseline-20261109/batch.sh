#!/bin/bash
#SBATCH --job-name=sharkmer-v4-baseline
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=tmp/v4-performance-baseline-20261109/logs/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module reset
source tmp/uncertain-products-20260927/toolchain/activate.sh
export CARGO_TARGET_DIR="$PWD/tmp/v4-performance-baseline-20261109/target"
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=1 LC_ALL=C
export PATH="/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin:$PATH"
export PYTHONPATH="$PWD/tmp/prerelease-checks-20261010/python-deps${PYTHONPATH:+:$PYTHONPATH}"
root="$PWD"
cargo build --manifest-path tmp/v4-performance-baseline-20261109/source/Cargo.toml --release --locked --offline --bin sharkmer
mkdir -p tmp/v4-performance-baseline-20261109/bin
cp "$CARGO_TARGET_DIR/release/sharkmer" tmp/v4-performance-baseline-20261109/bin/sharkmer-v3.2.1
sha256sum tmp/v4-performance-baseline-20261109/bin/sharkmer-v3.2.1
python benchmarks/benchmark_results/v4-performance-baseline-20261109/run_baseline.py --binary "$root/tmp/v4-performance-baseline-20261109/bin/sharkmer-v3.2.1" --output "$root/tmp/v4-performance-baseline-20261109/run-v1"
