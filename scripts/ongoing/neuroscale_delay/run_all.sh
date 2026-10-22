#!/bin/bash
# Reproduce the actsim delay analysis end to end.
#
# Usage: scripts/ongoing/neuroscale_delay/run_all.sh <actsim_ticks_core_timing.csv> [out_dir]
#
# Needs: the SANA-FE Python module importable (pip install . or PYTHONPATH to a
# build dir), the network JSON at runs/neuroscale/1.json (or NEUROSCALE_JSON),
# and pandas, numpy, matplotlib, pyyaml. Runs from the project root.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
ACTSIM_CSV="$(realpath "$1")"
OUT="$(realpath -m "${2:-$ROOT/runs/neuroscale/delay_analysis}")"
DATA="$OUT/data"; FIGS="$OUT/figs"
mkdir -p "$DATA" "$FIGS"
export MPLBACKEND=Agg NEUROSCALE_LOG_DIR="$DATA"
cd "$ROOT"
run() { echo; echo "===== $1"; python3 "$HERE/$1" "${@:2}" 2>&1 | grep -v '^\[' ; }
run build_dataset.py "$ROOT" "${NEUROSCALE_ARCH:-$ROOT/arch/neuroscale3.yaml}" "$ACTSIM_CSV" "$DATA"
run decompose.py "$DATA" 6
run correlate.py "$DATA"
run predictors.py "$DATA"
run stalls.py "$DATA" 3
run tile_level.py "$DATA" "$FIGS"
run receiver.py "$DATA" "$FIGS"
run release.py "$DATA" "$FIGS"
run backlog.py "$DATA" "$FIGS"
run stall_rule.py "$DATA" "$FIGS"
run tile_sync.py "$DATA"
run forward_sim.py "$DATA" "$FIGS"
echo; echo "Outputs in $OUT"
