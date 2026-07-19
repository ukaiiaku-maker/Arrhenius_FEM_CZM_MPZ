#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
MANIFEST_ROOT="${MANIFEST_ROOT:-runs/mpz_v9_10_4_dbtt_quality_diversity_v1}"
OUTROOT="${OUTROOT:-runs/mpz_v9_10_4_dbtt_spatial_screen_50um_v1}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
MAX_PER_CLASS="${MAX_PER_CLASS:-10}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-50}"
DA_UM="${DA_UM:-5}"
DK="${DK:-0.5}"
KDOT="${KDOT:-0.005}"
KMAX="${KMAX:-80}"
MPZ_LENGTH_UM="${MPZ_LENGTH_UM:-50}"
MPZ_N_BINS="${MPZ_N_BINS:-80}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

# This runner changes only which v9.10.3 candidates are promoted. The spatial
# mechanics remains the validated v9.10.2 independent-shape isotropic MPZ.
conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
  promote_mpz_v9_10_2_spatial.py \
  --manifest-root "$MANIFEST_ROOT" \
  --classes DBTT \
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
