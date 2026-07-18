#!/usr/bin/env bash
set -euo pipefail

# v10.0.5.4 VHCF first-passage remote stress-range campaign.
#
# The local hazard/state increments determine cycle-block size. MAX_BLOCK_CYCLES
# defaults to infinity so the physical cycle horizon, not an artificial numerical
# ceiling, limits rare-event jumps after the transient process-zone state settles.

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/${CONDA_ENV}/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || command -v python)"
fi

MODE="${MODE:-smoke}"                    # smoke | pilot | full
MATERIAL="${MATERIAL:-DBTT}"
TEMPERATURES="${TEMPERATURES:-700}"
DELTA_SIGMA_MPA="${DELTA_SIGMA_MPA:-250 300 350}"
R="${R:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
OUTROOT="${OUTROOT:-runs/v10_0_5_4_${MATERIAL}_vhcf_${MODE}}"

case "$MODE" in
  smoke)
    CYCLES_MAX="${CYCLES_MAX:-1e5}"
    MAX_BLOCKS="${MAX_BLOCKS:-500}"
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
    CYCLES_MAX="${CYCLES_MAX:-1e14}"
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

MAX_BLOCK_CYCLES="${MAX_BLOCK_CYCLES:-inf}"
BLOCK_CYCLES="${BLOCK_CYCLES:-1e4}"
MIN_BLOCK_CYCLES="${MIN_BLOCK_CYCLES:-1e-6}"
TARGET_DB="${TARGET_DB:-0.01}"
TARGET_DN_STORE="${TARGET_DN_STORE:-0.01}"
TARGET_DN_EMIT="${TARGET_DN_EMIT:-0.10}"
TARGET_DN_MOBILE="${TARGET_DN_MOBILE:-0.10}"
TARGET_DN_ESCAPE="${TARGET_DN_ESCAPE:-0.10}"
TARGET_DN_PEIERLS="${TARGET_DN_PEIERLS:-inf}"
TARGET_DN_TAYLOR="${TARGET_DN_TAYLOR:-inf}"
N_PHASE="${N_PHASE:-96}"
CYCLIC_MECHANICS_PHASES="${CYCLIC_MECHANICS_PHASES:-16}"
RESOLVE_CYCLIC_MECHANICS="${RESOLVE_CYCLIC_MECHANICS:-0}"
FAIL_ON_CENSOR="${FAIL_ON_CENSOR:-0}"
PRINT_EVERY="${PRINT_EVERY:-10}"

if [[ -e "$OUTROOT" && "${KEEP_EXISTING:-0}" != "1" ]]; then
  echo "ERROR: output path already exists: $OUTROOT" >&2
  echo "Use a new versioned OUTROOT or set KEEP_EXISTING=1." >&2
  exit 2
fi

"$PYTHON_BIN" -m compileall -q arrhenius_fracture \
  run_v10_0_5_4_vhcf_delta_sigma.py \
  run_v10_0_5_4_vhcf_delta_sigma_compat.py

KEEP_FLAG=""
if [[ "${KEEP_EXISTING:-0}" == "1" ]]; then
  KEEP_FLAG="--keep-existing"
fi

CYCLIC_FLAG=""
if [[ "$RESOLVE_CYCLIC_MECHANICS" == "1" ]]; then
  CYCLIC_FLAG="--resolve-cyclic-mechanics"
fi

CENSOR_FLAG=""
if [[ "$FAIL_ON_CENSOR" == "1" ]]; then
  CENSOR_FLAG="--fail-on-censor"
fi

# shellcheck disable=SC2086
"$PYTHON_BIN" run_v10_0_5_4_vhcf_delta_sigma_compat.py \
  --out "$OUTROOT" \
  --material-class "$MATERIAL" \
  --temperatures $TEMPERATURES \
  --delta-sigma-MPa $DELTA_SIGMA_MPA \
  --R "$R" \
  --frequency-Hz "$FREQUENCY_HZ" \
  --cycles-max "$CYCLES_MAX" \
  --block-cycles "$BLOCK_CYCLES" \
  --max-block-cycles "$MAX_BLOCK_CYCLES" \
  --min-block-cycles "$MIN_BLOCK_CYCLES" \
  --max-blocks "$MAX_BLOCKS" \
  --target-dB "$TARGET_DB" \
  --target-dN-store "$TARGET_DN_STORE" \
  --target-dN-emit "$TARGET_DN_EMIT" \
  --target-dN-mobile "$TARGET_DN_MOBILE" \
  --target-dN-escape "$TARGET_DN_ESCAPE" \
  --target-dN-peierls "$TARGET_DN_PEIERLS" \
  --target-dN-taylor "$TARGET_DN_TAYLOR" \
  --n-phase "$N_PHASE" \
  --cyclic-mechanics-phases "$CYCLIC_MECHANICS_PHASES" \
  --target-extension-um "$TARGET_EXTENSION_UM" \
  --nx "$NX" --ny "$NY" \
  --save-snapshots "$SAVE_SNAPSHOTS" \
  --print-every "$PRINT_EVERY" \
  ${CYCLIC_FLAG:+$CYCLIC_FLAG} \
  ${CENSOR_FLAG:+$CENSOR_FLAG} \
  ${KEEP_FLAG:+$KEEP_FLAG}

cat <<EOF
v10.0.5.4 VHCF fatigue campaign complete
out=$OUTROOT
material=$MATERIAL
T=[$TEMPERATURES]
DeltaSigma_MPa=[$DELTA_SIGMA_MPA]
R=$R
frequency_Hz=$FREQUENCY_HZ
cycles_max=$CYCLES_MAX
max_block_cycles=$MAX_BLOCK_CYCLES
resolve_cyclic_mechanics=$RESOLVE_CYCLIC_MECHANICS
summary=$OUTROOT/K_vs_delta_sigma.csv
block_diagnostics=$OUTROOT/fatigue_block_diagnostics_v10_0_5_4.csv
EOF
