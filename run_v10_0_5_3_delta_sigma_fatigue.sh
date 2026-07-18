#!/usr/bin/env bash
set -euo pipefail

# v10.0.5.3 remote stress-range fatigue campaign.
# The runner calibrates displacement from the model reaction force, then reports
# KJmax/DeltaKJ versus Delta-sigma and da/dN versus DeltaKJ.

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/${CONDA_ENV}/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || command -v python)"
fi

MODE="${MODE:-smoke}"                    # smoke | pilot | full
MATERIAL="${MATERIAL:-DBTT}"
TEMPERATURES="${TEMPERATURES:-300 500 700 900}"
DELTA_SIGMA_MPA="${DELTA_SIGMA_MPA:-200 250 300 350 400}"
R="${R:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
OUTROOT="${OUTROOT:-runs/v10_0_5_3_${MATERIAL}_delta_sigma_fatigue_${MODE}}"

case "$MODE" in
  smoke)
    CYCLES_MAX="${CYCLES_MAX:-1e5}"
    MAX_BLOCKS="${MAX_BLOCKS:-200}"
    TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-5}"
    NX="${NX:-40}"; NY="${NY:-80}"
    SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-0}"
    ;;
  pilot)
    CYCLES_MAX="${CYCLES_MAX:-1e9}"
    MAX_BLOCKS="${MAX_BLOCKS:-4000}"
    TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-50}"
    NX="${NX:-60}"; NY="${NY:-120}"
    SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-6}"
    ;;
  full)
    CYCLES_MAX="${CYCLES_MAX:-1e12}"
    MAX_BLOCKS="${MAX_BLOCKS:-20000}"
    TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-250}"
    NX="${NX:-80}"; NY="${NY:-160}"
    SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-12}"
    ;;
  *)
    echo "ERROR: MODE must be smoke, pilot, or full" >&2
    exit 2
    ;;
esac

if [[ -e "$OUTROOT" && "${KEEP_EXISTING:-0}" != "1" ]]; then
  echo "ERROR: output path already exists: $OUTROOT" >&2
  echo "Use a new versioned OUTROOT or set KEEP_EXISTING=1." >&2
  exit 2
fi

"$PYTHON_BIN" -m compileall -q arrhenius_fracture \
  run_v10_0_5_3_delta_sigma_fatigue.py

# Use a scalar rather than an empty array for compatibility with the macOS
# Bash 3.2 nounset behavior under `set -u`.
KEEP_FLAG=""
if [[ "${KEEP_EXISTING:-0}" == "1" ]]; then
  KEEP_FLAG="--keep-existing"
fi

# shellcheck disable=SC2086
"$PYTHON_BIN" run_v10_0_5_3_delta_sigma_fatigue.py \
  --out "$OUTROOT" \
  --material-class "$MATERIAL" \
  --temperatures $TEMPERATURES \
  --delta-sigma-MPa $DELTA_SIGMA_MPA \
  --R "$R" \
  --frequency-Hz "$FREQUENCY_HZ" \
  --cycles-max "$CYCLES_MAX" \
  --max-blocks "$MAX_BLOCKS" \
  --target-extension-um "$TARGET_EXTENSION_UM" \
  --nx "$NX" --ny "$NY" \
  --save-snapshots "$SAVE_SNAPSHOTS" \
  ${KEEP_FLAG:+$KEEP_FLAG}

cat <<EOF
v10.0.5.3 fatigue campaign complete
out=$OUTROOT
material=$MATERIAL
T=[$TEMPERATURES]
DeltaSigma_MPa=[$DELTA_SIGMA_MPA]
R=$R
frequency_Hz=$FREQUENCY_HZ
summary=$OUTROOT/K_vs_delta_sigma.csv
EOF
