#!/bin/bash
set -euo pipefail

# $SRC_DIR is the tarball root — Cargo.toml is at the top level
cd "${SRC_DIR}"

cargo install --locked --no-track --root "${PREFIX}" --path .
