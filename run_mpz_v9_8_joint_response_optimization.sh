#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTROOT="${OUTROOT:-runs/mpz_v9_8_joint_response_optimization_v1}"
TARGET_CLASSES="${TARGET_CLASSES:-ceramic weakT DBTT}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
SEED_COUNT="${SEED_COUNT:-6}"
SEED_POOL_SIZE="${SEED_POOL_SIZE:-120}"
DE_MAXITER="${DE_MAXITER:-40}"
DE_POPSIZE="${DE_POPSIZE:-8}"
LOCAL_MAXITER="${LOCAL_MAXITER:-300}"
DK="${DK:-0.5}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-65}"
SHORTLIST_COUNT="${SHORTLIST_COUNT:-20}"
MAX_JOBS="${MAX_JOBS:-2}"
SEED="${SEED:-98017}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

run_one() {
  local class_name="$1"
  local class_seed="$2"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting class=$class_name"
  conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
    optimize_mpz_v9_8_joint_response.py \
    --target-class "$class_name" \
    --temperatures "$TEMPERATURES" \
    --seed-count "$SEED_COUNT" \
    --seed-pool-size "$SEED_POOL_SIZE" \
    --de-maxiter "$DE_MAXITER" \
    --de-popsize "$DE_POPSIZE" \
    --local-maxiter "$LOCAL_MAXITER" \
    --dK "$DK" \
    --Kdot "$KDOT" \
    --Kmax "$KMAX" \
    --shortlist-count "$SHORTLIST_COUNT" \
    --seed "$class_seed" \
    --resume \
    --out "$OUTROOT"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] completed class=$class_name"
}

pids=()
index=0
for class_name in $TARGET_CLASSES; do
  while [ "${#pids[@]}" -ge "$MAX_JOBS" ]; do
    new=()
    for pid in "${pids[@]}"; do
      if kill -0 "$pid" 2>/dev/null; then
        new+=("$pid")
      else
        wait "$pid"
      fi
    done
    pids=("${new[@]}")
    [ "${#pids[@]}" -lt "$MAX_JOBS" ] || sleep 2
  done
  run_one "$class_name" "$((SEED + 10000 * index))" &
  pids+=("$!")
  index=$((index + 1))
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] all v9.8 joint optimizations complete"
