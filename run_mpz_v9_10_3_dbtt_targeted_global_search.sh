#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTROOT="${OUTROOT:-runs/mpz_v9_10_3_dbtt_targeted_global_search_v1}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
RESTARTS="${RESTARTS:-3}"
DE_MAXITER="${DE_MAXITER:-45}"
DE_POPSIZE="${DE_POPSIZE:-6}"
LOCAL_MAXITER="${LOCAL_MAXITER:-250}"
DK="${DK:-0.5}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-80}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-500}"
DA_UM="${DA_UM:-5}"
SEED="${SEED:-910317}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting v9.10.3 target-aware DBTT search"
conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
  optimize_mpz_v9_10_3_dbtt_targeted_global.py \
  --target-class DBTT \
  --temperatures "$TEMPERATURES" \
  --restarts "$RESTARTS" \
  --de-maxiter "$DE_MAXITER" \
  --de-popsize "$DE_POPSIZE" \
  --local-maxiter "$LOCAL_MAXITER" \
  --seed "$SEED" \
  --dK "$DK" \
  --Kdot "$KDOT" \
  --Kmax "$KMAX" \
  --target-extension-um "$TARGET_EXTENSION_UM" \
  --da-um "$DA_UM" \
  --out "$OUTROOT"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] completed v9.10.3 target-aware DBTT search"
