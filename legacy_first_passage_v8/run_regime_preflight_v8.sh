#!/usr/bin/env bash
set -euo pipefail
# Fastest useful check: all four regimes at 300 K for one orientation.
# Usage: THETA=30 STEPS=3000 bash run_regime_preflight_v8.sh preflight_theta30
OUT=${1:-regime_preflight_v8_theta${THETA:-30}}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TEMPS="300" REGIMES="ceramic peak weakT DBTT" STEPS=${STEPS:-3000} THETA=${THETA:-30} \
  bash "$SCRIPT_DIR/run_regime_sweep_v8_single_orientation.sh" "$OUT"
