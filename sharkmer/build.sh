#!/bin/bash
set -euo pipefail

# $SRC_DIR is the tarball root (sharkmer-1.0.0/)
# The Rust crate is one level deeper
cd "${SRC_DIR}/sharkmer"

cargo-bundle-licenses --format json --output licenses/LICENSE.dependencies.json

cargo install --locked --no-track --root "${PREFIX}" --path .
