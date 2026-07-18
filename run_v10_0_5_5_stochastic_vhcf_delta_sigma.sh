#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/${CONDA_ENV}/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || command -v python)"
fi

MODE="${MODE:-smoke}"                    # smoke | pilot | full
MATERIAL="${MATERIAL:-DBTT}"
TEMPERATURES="${TEMPERATURES:-700}"
DELTA_SIGMA_MPA="${DELTA_SIGMA_MPA:-350}"
R="${R:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
STOCHASTIC_SEED="${STOCHASTIC_SEED:-1}"
OUTROOT="${OUTROOT:-runs/v10_0_5_5_${MATERIAL}_stochastic_vhcf_${MODE}_seed${STOCHASTIC_SEED}}"

case "$MODE" in
  smoke)
    CYCLES_MAX="${CYCLES_MAX:-1e5}"
    MAX_BLOCKS="${MAX_BLOCKS:-1000}"
    TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-5}"
    NX="${NX:-40}"; NY="${NY:-80}"; SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-0}"
    ;;
  pilot)
    CYCLES_MAX="${CYCLES_MAX:-1e8}"
    MAX_BLOCKS="${MAX_BLOCKS:-5000}"
    TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-50}"
    NX="${NX:-60}"; NY="${NY:-120}"; SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-6}"
    ;;
  full)
    # 1e14 remains supported through CYCLES_MAX, but is not a mandatory target.
    CYCLES_MAX="${CYCLES_MAX:-1e12}"
    MAX_BLOCKS="${MAX_BLOCKS:-20000}"
    TARGET_EXTENSION_UM="${TARGET_EXTENSION_UM:-250}"
    NX="${NX:-80}"; NY="${NY:-160}"; SAVE_SNAPSHOTS="${SAVE_SNAPSHOTS:-12}"
    ;;
  *) echo "ERROR: MODE must be smoke, pilot, or full" >&2; exit 2 ;;
esac

MAX_BLOCK_CYCLES="${MAX_BLOCK_CYCLES:-inf}"
BLOCK_CYCLES="${BLOCK_CYCLES:-1e4}"
MIN_BLOCK_CYCLES="${MIN_BLOCK_CYCLES:-1e-6}"
TARGET_DB="${TARGET_DB:-0.01}"
TARGET_DN_STORE="${TARGET_DN_STORE:-0.05}"
TARGET_DN_EMIT="${TARGET_DN_EMIT:-inf}"
TARGET_DN_MOBILE="${TARGET_DN_MOBILE:-inf}"
TARGET_DN_ESCAPE="${TARGET_DN_ESCAPE:-0.25}"
TARGET_DN_PEIERLS="${TARGET_DN_PEIERLS:-inf}"
TARGET_DN_TAYLOR="${TARGET_DN_TAYLOR:-inf}"
N_PHASE="${N_PHASE:-96}"
RESOLVE_CYCLIC_MECHANICS="${RESOLVE_CYCLIC_MECHANICS:-0}"
FAIL_ON_CENSOR="${FAIL_ON_CENSOR:-0}"
PRINT_EVERY="${PRINT_EVERY:-10}"

EVENT_STATISTICS="${EVENT_STATISTICS:-stochastic}"
STOCHASTIC_EMISSION="${STOCHASTIC_EMISSION:-1}"
STOCHASTIC_BLOCKS="${STOCHASTIC_BLOCKS:-1}"
RARE_EVENT_TARGET="${RARE_EVENT_TARGET:-0.25}"
TAU_LEAP_TARGET="${TAU_LEAP_TARGET:-3.0}"
TAU_SWITCH_EXPECTED_EVENTS="${TAU_SWITCH_EXPECTED_EVENTS:-10.0}"
# Cache remains opt-in until a cache-on/cache-off equivalence smoke passes.
VHCF_FEM_CACHE="${VHCF_FEM_CACHE:-0}"

if [[ -e "$OUTROOT" && "${KEEP_EXISTING:-0}" != "1" ]]; then
  echo "ERROR: output path already exists: $OUTROOT" >&2
  echo "Use a new versioned OUTROOT or set KEEP_EXISTING=1." >&2
  exit 2
fi

export ARRHENIUS_EVENT_STATISTICS="$EVENT_STATISTICS"
export ARRHENIUS_STOCHASTIC_EMISSION="$STOCHASTIC_EMISSION"
export ARRHENIUS_STOCHASTIC_SEED="$STOCHASTIC_SEED"
export ARRHENIUS_STOCHASTIC_BLOCKS="$STOCHASTIC_BLOCKS"
export ARRHENIUS_RARE_EVENT_TARGET="$RARE_EVENT_TARGET"
export ARRHENIUS_TAU_LEAP_TARGET="$TAU_LEAP_TARGET"
export ARRHENIUS_TAU_SWITCH_EXPECTED_EVENTS="$TAU_SWITCH_EXPECTED_EVENTS"
export ARRHENIUS_VHCF_FEM_CACHE="$VHCF_FEM_CACHE"

"$PYTHON_BIN" -m compileall -q arrhenius_fracture \
  run_v10_0_5_5_stochastic_vhcf_delta_sigma.py \
  run_v10_0_5_5_stochastic_vhcf_delta_sigma_compat.py

KEEP_FLAG=""; [[ "${KEEP_EXISTING:-0}" == "1" ]] && KEEP_FLAG="--keep-existing"
CYCLIC_FLAG=""; [[ "$RESOLVE_CYCLIC_MECHANICS" == "1" ]] && CYCLIC_FLAG="--resolve-cyclic-mechanics"
CENSOR_FLAG=""; [[ "$FAIL_ON_CENSOR" == "1" ]] && CENSOR_FLAG="--fail-on-censor"

# shellcheck disable=SC2086
"$PYTHON_BIN" run_v10_0_5_5_stochastic_vhcf_delta_sigma_compat.py \
  --out "$OUTROOT" \
  --material-class "$MATERIAL" \
  --temperatures $TEMPERATURES \
  --delta-sigma-MPa $DELTA_SIGMA_MPA \
  --R "$R" --frequency-Hz "$FREQUENCY_HZ" \
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
  --target-extension-um "$TARGET_EXTENSION_UM" \
  --nx "$NX" --ny "$NY" \
  --save-snapshots "$SAVE_SNAPSHOTS" \
  --print-every "$PRINT_EVERY" \
  ${CYCLIC_FLAG:+$CYCLIC_FLAG} \
  ${CENSOR_FLAG:+$CENSOR_FLAG} \
  ${KEEP_FLAG:+$KEEP_FLAG}

cat <<EOF
v10.0.5.5 stochastic VHCF campaign complete
out=$OUTROOT
cycles_max=$CYCLES_MAX
stochastic_seed=$STOCHASTIC_SEED
rare_event_target=$RARE_EVENT_TARGET
tau_leap_target=$TAU_LEAP_TARGET
fem_cache=$VHCF_FEM_CACHE
summary=$OUTROOT/K_vs_delta_sigma.csv
EOF
