#!/usr/bin/env bash
# Validate all panel YAML files in panels/ against their declared samples.
# Requires the sharkmer-bench conda environment to be active.
#
# Usage:
#   conda activate sharkmer-bench
#   bash scripts/validate_all_panels.sh [--no-blast]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PANELS_DIR="$REPO_ROOT/panels"
SCRIPT="$REPO_ROOT/scripts/validate_panel.py"

# Pass through any arguments (e.g. --no-blast) to validate_panel.py.
EXTRA_ARGS=("$@")

panels=("$PANELS_DIR"/*.yaml)

if [ ${#panels[@]} -eq 0 ]; then
    echo "No panel YAML files found in $PANELS_DIR"
    exit 1
fi

echo "Validating ${#panels[@]} panels..."
echo ""

failed=()
succeeded=()

for panel in "${panels[@]}"; do
    name="$(basename "$panel")"
    echo "=== $name ==="
    if python "$SCRIPT" "$panel" "${EXTRA_ARGS[@]}"; then
        succeeded+=("$name")
    else
        failed+=("$name")
    fi
    echo ""
done

echo "==============================="
echo "Results: ${#succeeded[@]} passed, ${#failed[@]} failed out of ${#panels[@]} panels"
if [ ${#succeeded[@]} -gt 0 ]; then
    echo "  Passed: ${succeeded[*]}"
fi
if [ ${#failed[@]} -gt 0 ]; then
    echo "  Failed: ${failed[*]}"
    exit 1
fi
