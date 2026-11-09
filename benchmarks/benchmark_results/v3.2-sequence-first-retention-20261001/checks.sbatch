#!/bin/bash
#SBATCH --job-name=sequence-retention-checks
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=04:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
source tmp/uncertain-products-20260927/toolchain/activate.sh
export CARGO_TARGET_DIR="$PWD/tmp/sequence-retention-20261001/target"
mkdir -p tmp/sequence-retention-20261001/checks-v1
sha256sum Cargo.toml Cargo.lock src/pcr/mod.rs src/pcr/paths.rs src/pcr/uncertain.rs scripts/sharkmer_validate/uncertain.py tests/uncertain_products.rs tests/test_uncertain_output.py > tmp/sequence-retention-20261001/checks-v1/sources.before.sha256
cargo test --release --locked --offline > tmp/sequence-retention-20261001/checks-v1/test-default.log 2>&1
cargo test --release --locked --offline --no-default-features --features fxhashmap > tmp/sequence-retention-20261001/checks-v1/test-fxhash.log 2>&1
cargo clippy --all-targets --locked --offline -- -D warnings > tmp/sequence-retention-20261001/checks-v1/clippy-default.log 2>&1
cargo clippy --all-targets --locked --offline --no-default-features --features fxhashmap -- -D warnings > tmp/sequence-retention-20261001/checks-v1/clippy-fxhash.log 2>&1
cargo build --release --locked --offline > tmp/sequence-retention-20261001/checks-v1/build.log 2>&1
sha256sum Cargo.toml Cargo.lock src/pcr/mod.rs src/pcr/paths.rs src/pcr/uncertain.rs scripts/sharkmer_validate/uncertain.py tests/uncertain_products.rs tests/test_uncertain_output.py > tmp/sequence-retention-20261001/checks-v1/sources.after.sha256
cmp tmp/sequence-retention-20261001/checks-v1/sources.before.sha256 tmp/sequence-retention-20261001/checks-v1/sources.after.sha256
sha256sum "$CARGO_TARGET_DIR/release/sharkmer" > tmp/sequence-retention-20261001/checks-v1/binary.sha256
