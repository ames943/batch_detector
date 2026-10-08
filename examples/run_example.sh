#!/bin/bash
# Run batch_detector on the synthetic example data.
# Usage: bash examples/run_example.sh     (from anywhere)
# Set PYTHON to use a different interpreter, e.g. PYTHON=python3.11 bash examples/run_example.sh

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
PYTHON="${PYTHON:-python}"

# --n-grid leaves out n = 480 to keep the run short; everything else is default
"$PYTHON" "$ROOT/batch_detector.py" \
    --input  "$HERE/example_features.tsv" \
    --labels "$HERE/example_metadata.tsv" \
    --batch  cohort \
    --clr \
    --n-grid 60,120,240 \
    --output "$ROOT/outputs/example"

echo
echo "Done. Results are in $ROOT/outputs/example"
