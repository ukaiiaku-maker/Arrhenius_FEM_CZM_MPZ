#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
SELECTED_ROWS="${SELECTED_ROWS:-runs/mpz_v9_4_peierls_taylor_search_v1/strict_common_selection/pt_v9_4_recommended_intrinsic_rows.csv}"
TARGETS="${TARGETS:-mpz_three_class_design_targets.csv}"
OUTROOT="${OUTROOT:-runs/mpz_v9_4_developed_state_search_v1}"
SAMPLES="${SAMPLES:-128}"
MAX_WORKERS="${MAX_WORKERS:-2}"
TEMPERATURES="${TEMPERATURES:-ceramic:300,900,1200;weakT:300,700,1200;DBTT:300,700,900,1200}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-400}"
DA_UM="${DA_UM:-5}"
DK="${DK:-0.25}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-65}"
TOP_COUNT="${TOP_COUNT:-12}"
EVENT_TOP_COUNT="${EVENT_TOP_COUNT:-3}"
SEED="${SEED:-94131}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

stamp() { date '+%Y-%m-%d %H:%M:%S'; }
run_python() {
  if command -v conda >/dev/null 2>&1; then
    conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u "$@"
  else
    "$PYTHON_BIN" -u "$@"
  fi
}

echo "[$(stamp)] MPZ v9.4 common developed-state search"
echo "[$(stamp)] selected_rows=$SELECTED_ROWS"
echo "[$(stamp)] targets=$TARGETS"
echo "[$(stamp)] output=$OUTROOT"
echo "[$(stamp)] samples=$SAMPLES workers=$MAX_WORKERS"
echo "[$(stamp)] temperatures=$TEMPERATURES"
echo "[$(stamp)] extension_um=$TARGET_EXTENSION_UM da_um=$DA_UM dK=$DK Kdot=$KDOT Kmax=$KMAX"

run_python search_mpz_v9_4_developed_state.py \
  --selected-rows "$SELECTED_ROWS" \
  --targets "$TARGETS" \
  --samples "$SAMPLES" \
  --seed "$SEED" \
  --max-workers "$MAX_WORKERS" \
  --temperatures "$TEMPERATURES" \
  --target-extension-um "$TARGET_EXTENSION_UM" \
  --da-um "$DA_UM" \
  --dK "$DK" \
  --Kdot "$KDOT" \
  --Kmax "$KMAX" \
  --top-count "$TOP_COUNT" \
  --event-top-count "$EVENT_TOP_COUNT" \
  --resume \
  --out "$OUTROOT"

echo "[$(stamp)] MPZ v9.4 common developed-state search complete"
