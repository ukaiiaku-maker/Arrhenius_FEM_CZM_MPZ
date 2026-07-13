#!/usr/bin/env bash
set -euo pipefail

# Long-growth 2-D material-response runs for the v8 sharp-front fatigue model.
# This is the second-stage atlas: it reruns only the growing/near-cutoff points,
# lets cracks grow to a target extension, and records snapshots by crack extension.
#
# Useful environment overrides:
#   OUTROOT=runs/v8_material_response_long_growth_2d
#   CASE_FILTER="FCC_like_case29 plastic_shielded_case64_M1"
#   KLIST_OVERRIDE="7.0 8.0"
#   TARGET_EXT_UM=150
#   BLOCKS=1500
#   CYCLES_MAX=1e12
#   MAKE_2D_PLOTS=0   # CSV-only/debug

python -m compileall -q arrhenius_fracture

OUTROOT="${OUTROOT:-runs/v8_material_response_long_growth_2d}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
CYCLES_MAX="${CYCLES_MAX:-1e12}"
BLOCKS="${BLOCKS:-1500}"
NX="${NX:-24}"
NY="${NY:-48}"
SNAPSHOTS="${SNAPSHOTS:-36}"
SNAPSHOT_COLS="${SNAPSHOT_COLS:-6}"
TARGET_EXT_UM="${TARGET_EXT_UM:-150}"
SNAPSHOT_BY_EXT_UM="${SNAPSHOT_BY_EXT_UM:-25}"
MAX_DA_PER_BLOCK_UM="${MAX_DA_PER_BLOCK_UM:-10}"
TARGET_DA_PER_BLOCK_UM="${TARGET_DA_PER_BLOCK_UM:-5}"
DA_PHYS="${DA_PHYS:-5e-6}"
TARGET_DB="${TARGET_DB:-0.02}"
TARGET_DN_STORE="${TARGET_DN_STORE:-0.01}"
TARGET_DN_EMIT="${TARGET_DN_EMIT:-0.20}"
TARGET_DN_MOBILE="${TARGET_DN_MOBILE:-0.20}"
CYCLE_PHASES="${CYCLE_PHASES:-8}"
CASE_FILTER="${CASE_FILTER:-}"
KLIST_OVERRIDE="${KLIST_OVERRIDE:-}"
MAKE_2D_PLOTS="${MAKE_2D_PLOTS:-1}"

should_run_case() {
  local label="$1"
  if [[ -z "${CASE_FILTER}" ]]; then
    return 0
  fi
  # CASE_FILTER may be comma- or space-separated.
  local normalized=" ${CASE_FILTER//,/ } "
  [[ "${normalized}" == *" ${label} "* ]]
}

run_case() {
  local label="$1"
  local default_klist="$2"
  local G00="$3"
  local sigc="$4"
  local expa="$5"
  local expn="$6"
  local floor="$7"
  local emitS="$8"
  local peierlsS="$9"
  local taylorS="${10}"

  if ! should_run_case "${label}"; then
    echo "=== Skipping ${label} due to CASE_FILTER ==="
    return 0
  fi

  local klist="${KLIST_OVERRIDE:-${default_klist}}"
  local plot_flag="--make-2d-plots"
  if [[ "${MAKE_2D_PLOTS}" == "0" || "${MAKE_2D_PLOTS}" == "false" ]]; then
    plot_flag="--no-make-2d-plots"
  fi

  echo "=== Long-growth 2-D run: ${label}; Kmax=[${klist}] ==="
  python run_v8_compare_1d_2d_K_sweep.py \
    --out "${OUTROOT}/${label}" \
    --Kmax-MPa-sqrt-m ${klist} \
    --T "${T}" \
    --R "${R}" \
    --frequency-Hz "${FREQ}" \
    --blocks "${BLOCKS}" \
    --cycles-max "${CYCLES_MAX}" \
    --block-cycles 1e5 \
    --max-block-cycles inf \
    --cycle-block-mode hazard_limited \
    --target-dB "${TARGET_DB}" \
    --target-dN-store "${TARGET_DN_STORE}" \
    --target-dN-emit "${TARGET_DN_EMIT}" \
    --target-dN-mobile "${TARGET_DN_MOBILE}" \
    --storage-model escape_limited \
    --calibrate-2d-K \
    --K-calib-iters 3 \
    --K-calib-tol 5e-3 \
    --no-stop-after-first-2d-fire \
    --cyclic-mechanics-phases "${CYCLE_PHASES}" \
    --nx "${NX}" --ny "${NY}" \
    --tip-h-fine 1e-6 \
    --tip-ratio 1.25 \
    --da-phys "${DA_PHYS}" \
    --target-da-per-block-um "${TARGET_DA_PER_BLOCK_UM}" \
    --target-crack-extension-um "${TARGET_EXT_UM}" \
    --snapshot-by-crack-extension-um "${SNAPSHOT_BY_EXT_UM}" \
    --max-da-per-block-um "${MAX_DA_PER_BLOCK_UM}" \
    --save-snapshots "${SNAPSHOTS}" \
    --snapshot-cols "${SNAPSHOT_COLS}" \
    ${plot_flag} \
    --min-global-forward 0.05 \
    --cleave-barrier-kind exp_floor \
    --cleave-exp-T-mode mu_scale \
    --cleave-G00-eV "${G00}" \
    --cleave-sigc0-GPa "${sigc}" \
    --cleave-exp-a "${expa}" \
    --cleave-exp-n "${expn}" \
    --cleave-floor-frac "${floor}" \
    --emit-energy-scale 0.75 \
    --emit-entropy-scale "${emitS}" \
    --peierls-energy-scale 0.00375 \
    --peierls-entropy-scale "${peierlsS}" \
    --peierls-stress-scale 1.0 \
    --taylor-energy-scale 0.015 \
    --taylor-entropy-scale "${taylorS}" \
    --taylor-stress-scale 1.0
}

# Default K lists focus on previously measured growth windows plus near-cutoff points.
run_case "FCC_like_case29"             "6.0 7.0 8.0"       1.0 2.5 0.70 0.6 0.020 0.375 0.001875 0.0075
run_case "shifted_ductile_case64"      "6.5 7.5 8.0"       1.0 3.0 0.70 0.6 0.010 0.375 0.001875 0.0075
run_case "steep_cleavage_case35"       "7.0 8.0"           1.0 2.5 0.70 1.0 0.020 0.0   0.0      0.0
run_case "slow_threshold_case101"      "7.5 8.0"           1.0 3.5 0.70 0.6 0.020 0.0   0.0      0.0
run_case "higher_barrier_case171"      "7.0 8.0"           1.1 2.5 0.70 0.6 0.005 0.0   0.0      0.0
run_case "plastic_shielded_case64_M1"  "6.5 7.0 7.5 8.0"   1.0 3.0 0.70 0.6 0.010 0.75  0.00375  0.015

python postprocess_v8_material_response_atlas.py \
  --root "${OUTROOT}" \
  --case-table selected_2d_material_response_cases.csv \
  --R "${R}" \
  --cycles-max "${CYCLES_MAX}"
