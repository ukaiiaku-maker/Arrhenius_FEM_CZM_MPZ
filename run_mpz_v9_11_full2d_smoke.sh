#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
CLASS=${CLASS:-DBTT}
T_K=${T_K:-700}
OUTROOT=${OUTROOT:-runs/mpz_v9_11_full2d_${CLASS}_${T_K}K_smoke_v1}
RUN_SOLVER=${RUN_SOLVER:-1}
NX=${NX:-24}
NY=${NY:-48}
STEPS=${STEPS:-1200}
DU=${DU:-2e-7}
DT=${DT:-8.4}
PRINT_EVERY=${PRINT_EVERY:-50}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-0}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-45}

class_key=$(printf '%s' "$CLASS" | tr '[:upper:]' '[:lower:]')
case "$class_key" in
  ceramic)
    CLASS_CANON=ceramic
    ;;
  weakt|weak_t|weak)
    CLASS_CANON=weakT
    ;;
  dbtt)
    CLASS_CANON=DBTT
    ;;
  *)
    echo "ERROR: CLASS must be ceramic, weakT, or DBTT; got '$CLASS'" >&2
    exit 2
    ;;
esac

MANIFEST="$PARAMETER_ROOT/$CLASS_CANON/spatial_promotion_manifest.csv"
if [[ ! -s "$MANIFEST" ]]; then
  echo "ERROR: selected material manifest not found: $MANIFEST" >&2
  exit 2
fi

python verify_mpz_v9_11_install.py .
python verify_mpz_v9_11_physics.py --parameter-root "$PARAMETER_ROOT"
python -m pytest -q \
  tests/test_mode_i_first_passage_v9_11.py \
  tests/test_mpz_v9_11_2d_coupling.py \
  tests/test_mpz_v9_10_2_independent_shapes.py \
  tests/test_bulk_pt_plasticity.py

if [[ "$RUN_SOLVER" != "1" ]]; then
  echo "Preflight passed; RUN_SOLVER=$RUN_SOLVER so the FEM/CZM solve was skipped."
  exit 0
fi

CASE_OUT="$OUTROOT/$CLASS_CANON"
mkdir -p "$CASE_OUT"

python -m arrhenius_fracture.mode_i_first_passage_v9_11 \
  --mpz-material-manifest "$MANIFEST" \
  --mpz-material-class "$CLASS_CANON" \
  --mpz-length-um "$MPZ_LENGTH_UM" \
  --mpz-n-bins "$MPZ_N_BINS" \
  --mpz-profile-sector-half-angle-deg 45 \
  --mpz-profile-damage-cutoff 0.85 \
  --mode 2d \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine 3e-6 --tip-ratio 1.25 \
  --dU "$DU" --dt "$DT" --steps "$STEPS" \
  --n-stagger 2 --print-every "$PRINT_EVERY" \
  --stop-after-first-fire --max-fronts 1 \
  --adaptive-events --adaptive-event-target .25 \
  --adaptive-min-frac 1e-8 --adaptive-grow 4 \
  --da-phys 5e-6 \
  --j-decomposition cluster \
  --rJ-cluster 20e-6 \
  --rJ-outer 25e-6 \
  --temperatures "$T_K" \
  --crack-backend adaptive_czm \
  --czm-max-angle-error-deg 35 \
  --crystal-aniso --crystal-compete \
  --crystal-theta-deg "$CRYSTAL_THETA_DEG" \
  --crystal-C11 523e9 \
  --crystal-C12 203e9 \
  --crystal-C44 160e9 \
  --cleave-gamma-aniso 0.3 \
  --crystal-material w \
  --multihit-m 3 --multihit-tau 1e-6 \
  --sigma-cap-GPa 0 \
  --save-snapshots "$SAVE_SNAPSHOTS" \
  --no-plots \
  --out "$CASE_OUT"

echo "Mode-I v9.11 full 2-D smoke complete: $CASE_OUT"
