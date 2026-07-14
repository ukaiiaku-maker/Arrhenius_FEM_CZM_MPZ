#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
INPUT_ROOT="${INPUT_ROOT:-runs/mpz_v9_8_1_joint_response_optimization_v1}"
OUTROOT="${OUTROOT:-runs/mpz_v9_9_1_barrier_continuation_v1}"
CLASSES="${CLASSES:-ceramic weakT DBTT}"
SCALES="${SCALES:-1.0 0.8 0.6 0.4 0.3}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
CANDIDATES_PER_CLASS="${CANDIDATES_PER_CLASS:-3}"
LOCAL_MAXITER="${LOCAL_MAXITER:-400}"
DK="${DK:-0.5}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-65}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
  continue_mpz_v9_9_1_barrier_scale.py \
  --input-root "$INPUT_ROOT" \
  --classes "$CLASSES" \
  --scales "$SCALES" \
  --temperatures "$TEMPERATURES" \
  --candidates-per-class "$CANDIDATES_PER_CLASS" \
  --local-maxiter "$LOCAL_MAXITER" \
  --dK "$DK" \
  --Kdot "$KDOT" \
  --Kmax "$KMAX" \
  --out "$OUTROOT"
