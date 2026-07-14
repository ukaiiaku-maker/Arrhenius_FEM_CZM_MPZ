#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTROOT="${OUTROOT:-runs/mpz_v9_6_broad_dbtt_map_v1}"
PT_SAMPLES="${PT_SAMPLES:-64}"
MAX_INTRINSIC="${MAX_INTRINSIC:-0}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
RECOVERY_RATE_S="${RECOVERY_RATE_S:-1e-5}"
BACKSTRESS_UNIT="${BACKSTRESS_UNIT:-1.4}"
SEED="${SEED:-96061}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

run_python() {
  if command -v conda >/dev/null 2>&1; then
    conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u "$@"
  else
    "$PYTHON_BIN" -u "$@"
  fi
}

run_python search_mpz_v9_6_broad_dbtt_map.py \
  --pt-samples "$PT_SAMPLES" \
  --max-intrinsic "$MAX_INTRINSIC" \
  --temperatures "$TEMPERATURES" \
  --recovery-rate-s "$RECOVERY_RATE_S" \
  --backstress-unit-MPa-sqrt-m-per-sqrt-N "$BACKSTRESS_UNIT" \
  --seed "$SEED" \
  --out "$OUTROOT"
