#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTROOT="${OUTROOT:-runs/mpz_v9_10_1_shared_shape_global_search_v1}"
TARGET_CLASSES="${TARGET_CLASSES:-ceramic weakT DBTT}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
RESTARTS="${RESTARTS:-3}"
DE_MAXITER="${DE_MAXITER:-60}"
DE_POPSIZE="${DE_POPSIZE:-8}"
LOCAL_MAXITER="${LOCAL_MAXITER:-250}"
MAX_JOBS="${MAX_JOBS:-2}"
DK="${DK:-0.5}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-80}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-500}"
DA_UM="${DA_UM:-5}"
SEED="${SEED:-9101017}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

run_one() {
  local class_name="$1"
  local class_seed="$2"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting v9.10.1 shared-shape class=$class_name"
  conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
    optimize_mpz_v9_10_1_shared_shape_global.py \
    --target-class "$class_name" \
    --temperatures "$TEMPERATURES" \
    --restarts "$RESTARTS" \
    --de-maxiter "$DE_MAXITER" \
    --de-popsize "$DE_POPSIZE" \
    --local-maxiter "$LOCAL_MAXITER" \
    --seed "$class_seed" \
    --dK "$DK" \
    --Kdot "$KDOT" \
    --Kmax "$KMAX" \
    --target-extension-um "$TARGET_EXTENSION_UM" \
    --da-um "$DA_UM" \
    --out "$OUTROOT"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] completed v9.10.1 shared-shape class=$class_name"
}

pids=()
index=0
for class_name in $TARGET_CLASSES; do
  while [ "${#pids[@]}" -ge "$MAX_JOBS" ]; do
    next=()
    for pid in "${pids[@]}"; do
      if kill -0 "$pid" 2>/dev/null; then
        next+=("$pid")
      else
        wait "$pid"
      fi
    done
    pids=("${next[@]}")
    [ "${#pids[@]}" -lt "$MAX_JOBS" ] || sleep 2
  done
  run_one "$class_name" "$((SEED + 10000 * index))" &
  pids+=("$!")
  index=$((index + 1))
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] all v9.10.1 shared-shape searches complete"
