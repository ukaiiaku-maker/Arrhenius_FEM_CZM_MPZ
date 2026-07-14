#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
MANIFEST="${MANIFEST:-runs/mpz_v9_9_barrier_continuation_v1/spatial_promotion_manifest.csv}"
OUTROOT="${OUTROOT:-runs/mpz_v9_9_spatial_promotion_v1}"
CLASSES="${CLASSES:-ceramic weakT DBTT}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
MAX_PER_CLASS="${MAX_PER_CLASS:-2}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-500}"
DA_UM="${DA_UM:-5}"
DK="${DK:-0.25}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-80}"
MPZ_LENGTH_UM="${MPZ_LENGTH_UM:-100}"
MPZ_N_BINS="${MPZ_N_BINS:-200}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
  promote_mpz_v9_9_spatial.py \
  --manifest "$MANIFEST" \
  --classes "$CLASSES" \
  --temperatures "$TEMPERATURES" \
  --max-per-class "$MAX_PER_CLASS" \
  --target-extension-um "$TARGET_EXTENSION_UM" \
  --da-um "$DA_UM" \
  --dK "$DK" \
  --Kdot "$KDOT" \
  --Kmax "$KMAX" \
  --mpz-length-um "$MPZ_LENGTH_UM" \
  --mpz-n-bins "$MPZ_N_BINS" \
  --out "$OUTROOT"
