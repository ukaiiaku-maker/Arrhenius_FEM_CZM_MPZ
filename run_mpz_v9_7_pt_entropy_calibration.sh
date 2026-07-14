#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTROOT="${OUTROOT:-runs/mpz_v9_7_pt_entropy_calibration_v1}"
CLASSES="${CLASSES:-ceramic weakT DBTT}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
STRAIN_RATES="${STRAIN_RATES:-1e-5 1e-3}"
RHO_VALUES="${RHO_VALUES:-5e12 1e14 1e16}"
ENERGY_RATIO_POINTS="${ENERGY_RATIO_POINTS:-17}"
MAGNITUDE_TOP_PER_CLASS="${MAGNITUDE_TOP_PER_CLASS:-16}"
ENTROPY_SAMPLES="${ENTROPY_SAMPLES:-256}"
TARGET_REFERENCE_STRESS_GPA="${TARGET_REFERENCE_STRESS_GPA:-2.0}"
SEED="${SEED:-97017}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

run_python() {
  if command -v conda >/dev/null 2>&1; then
    conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u "$@"
  else
    "$PYTHON_BIN" -u "$@"
  fi
}

run_python calibrate_mpz_v9_7_pt_entropy.py \
  --classes "$CLASSES" \
  --temperatures "$TEMPERATURES" \
  --strain-rates "$STRAIN_RATES" \
  --rho-values-m2 "$RHO_VALUES" \
  --energy-ratio-points "$ENERGY_RATIO_POINTS" \
  --magnitude-top-per-class "$MAGNITUDE_TOP_PER_CLASS" \
  --entropy-samples "$ENTROPY_SAMPLES" \
  --target-reference-stress-GPa "$TARGET_REFERENCE_STRESS_GPA" \
  --seed "$SEED" \
  --out "$OUTROOT"
