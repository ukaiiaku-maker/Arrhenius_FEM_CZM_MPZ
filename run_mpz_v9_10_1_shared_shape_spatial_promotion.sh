#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
MANIFEST_ROOT="${MANIFEST_ROOT:-runs/mpz_v9_10_1_shared_shape_global_search_v1}"
OUTROOT="${OUTROOT:-runs/mpz_v9_10_1_shared_shape_spatial_promotion_v1}"
CLASSES="${CLASSES:-ceramic weakT DBTT}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
MAX_PER_CLASS="${MAX_PER_CLASS:-3}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-500}"
DA_UM="${DA_UM:-5}"
DK="${DK:-0.25}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-80}"
MPZ_LENGTH_UM="${MPZ_LENGTH_UM:-100}"
MPZ_N_BINS="${MPZ_N_BINS:-200}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

# v9.10 spatial promotion already uses the emission alpha/n values for the
# emission, Peierls, and Taylor surfaces.  The v9.10.1 manifest also sets the
# cleavage alpha/n values equal to those same shared values, so this promotes
# the exact one-shape parameterization without a second adapter.
conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
  promote_mpz_v9_10_spatial.py \
  --manifest-root "$MANIFEST_ROOT" \
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
