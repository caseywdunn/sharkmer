#!/bin/bash
set -euo pipefail

cd sharkmer

cargo-bundle-licenses --format json --output licenses/LICENSE.dependencies.json

cargo install --locked --no-track --root "$PREFIX" --path .
