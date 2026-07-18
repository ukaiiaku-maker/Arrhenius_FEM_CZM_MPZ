#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/${CONDA_ENV}/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || command -v python)"
fi

ACTION="${ACTION:-all}"                 # audit | bracket | all
MATERIAL="${MATERIAL:-DBTT}"
TEMPERATURE_K="${TEMPERATURE_K:-700}"
R="${R:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
STOCHASTIC_SEED="${STOCHASTIC_SEED:-1}"
CYCLES_MAX="${CYCLES_MAX:-1e7}"
MAX_BLOCKS="${MAX_BLOCKS:-10000}"
TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-5}"
OUTROOT="${OUTROOT:-runs/v10_0_5_6_${MATERIAL}_${TEMPERATURE_K}K_kj_audit_bracket_seed${STOCHASTIC_SEED}_v1}"

CONTOUR_OUTER_UM="${CONTOUR_OUTER_UM:-80 100 140 180 240 300 360 400}"
TARGET_KJMAX="${TARGET_KJMAX:-2 4 6 8 10 12 16 20 24}"
AUDIT_DELTA_SIGMA_MPA="${AUDIT_DELTA_SIGMA_MPA:-100}"
PLATEAU_REL_TOL="${PLATEAU_REL_TOL:-0.10}"
PLATEAU_MIN_POINTS="${PLATEAU_MIN_POINTS:-3}"
BISECTION_REFINEMENTS="${BISECTION_REFINEMENTS:-3}"
NX="${NX:-40}"
NY="${NY:-80}"
TIP_H_FINE_M="${TIP_H_FINE_M:-2.5e-6}"
TIP_RATIO="${TIP_RATIO:-1.2}"
N_PHASE="${N_PHASE:-96}"
PRINT_EVERY="${PRINT_EVERY:-10}"

if [[ -e "$OUTROOT" && "${KEEP_EXISTING:-0}" != "1" ]]; then
  echo "ERROR: output path already exists: $OUTROOT" >&2
  echo "Use a new versioned OUTROOT or set KEEP_EXISTING=1." >&2
  exit 2
fi
mkdir -p "$OUTROOT"

"$PYTHON_BIN" -m compileall -q \
  arrhenius_fracture/kj_audit_v10056.py \
  run_v10_0_5_6_stochastic_delta_sigma.py \
  run_v10_0_5_6_stochastic_delta_sigma_audited.py \
  run_v10_0_5_6_kj_audit_bracket.py \
  run_v10_0_5_6_kj_audit_bracket_audited.py

CMD=(
  "$PYTHON_BIN" run_v10_0_5_6_kj_audit_bracket_audited.py "$ACTION"
  --out "$OUTROOT"
  --material-class "$MATERIAL"
  --temperature-K "$TEMPERATURE_K"
  --R "$R"
  --frequency-Hz "$FREQUENCY_HZ"
  --stochastic-seed "$STOCHASTIC_SEED"
  --cycles-max "$CYCLES_MAX"
  --max-blocks "$MAX_BLOCKS"
  --target-extension-um "$TARGET_EXTENSION_UM"
  --audit-delta-sigma-MPa "$AUDIT_DELTA_SIGMA_MPA"
  --plateau-relative-tolerance "$PLATEAU_REL_TOL"
  --plateau-minimum-points "$PLATEAU_MIN_POINTS"
  --bisection-refinements "$BISECTION_REFINEMENTS"
  --nx "$NX" --ny "$NY"
  --tip-h-fine-m "$TIP_H_FINE_M"
  --tip-ratio "$TIP_RATIO"
  --n-phase "$N_PHASE"
  --print-every "$PRINT_EVERY"
  --contour-outer-um $CONTOUR_OUTER_UM
  --target-KJmax-MPa-sqrt-m $TARGET_KJMAX
)

if [[ "$ACTION" == "bracket" ]]; then
  SELECTED_CONTOUR_JSON="${SELECTED_CONTOUR_JSON:?Set SELECTED_CONTOUR_JSON for ACTION=bracket}"
  CMD+=(--selected-contour-json "$SELECTED_CONTOUR_JSON")
fi

"${CMD[@]}"

cat <<EOF
v10.0.5.6 KJ audit / first-passage bracket complete
out=$OUTROOT
action=$ACTION
selected_contour=$OUTROOT/selected_KJ_contour_v10_0_5_6.json
KJ_audit=$OUTROOT/KJ_contour_sweep_v10_0_5_6.csv
bracket=$OUTROOT/first_passage_bracket_v10_0_5_6.json
EOF
