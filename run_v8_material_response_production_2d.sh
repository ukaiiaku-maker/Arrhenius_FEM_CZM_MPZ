#!/usr/bin/env bash
set -euo pipefail

# Production 2-D material-response atlas for the v8 sharp-front fatigue model.
#
# Goal:
#   Capture the full experimental crack-growth envelope for the selected material
#   response classes: low-DeltaK no measurable growth, Paris-like growth, strong
#   shielding/cutoff behavior, and high-DeltaK rapid fracture (<1 cycle when it occurs).
#
# This is intentionally broader than run_v8_material_response_long_growth_2d.sh.
# It sweeps low to very high Kmax and uses adaptive cycle stepping with a large
# physical cycle horizon. High-K points finish quickly by target crack extension;
# low-K points become censored upper bounds only after the cycle horizon/block budget.
#
# Useful environment overrides:
#   OUTROOT=runs/v8_material_response_production_2d
#   CASE_FILTER="FCC_like_case29 plastic_shielded_case64_M1"
#   KLIST_OVERRIDE="4.0 5.0 6.0 7.0 8.0 10.0 12.0"
#   PRODUCTION_LEVEL=coarse|full
#   TARGET_EXT_UM=250
#   CYCLES_MAX=1e12
#   BLOCKS=3000
#   MAKE_2D_PLOTS=0     # CSV-only/debug
#   NX=24 NY=48

python -m compileall -q arrhenius_fracture

OUTROOT="${OUTROOT:-runs/v8_material_response_production_2d}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
CYCLES_MAX="${CYCLES_MAX:-1e12}"
BLOCKS="${BLOCKS:-3000}"
NX="${NX:-24}"
NY="${NY:-48}"
SNAPSHOTS="${SNAPSHOTS:-48}"
SNAPSHOT_COLS="${SNAPSHOT_COLS:-6}"
TARGET_EXT_UM="${TARGET_EXT_UM:-250}"
SNAPSHOT_BY_EXT_UM="${SNAPSHOT_BY_EXT_UM:-25}"
MAX_DA_PER_BLOCK_UM="${MAX_DA_PER_BLOCK_UM:-10}"
TARGET_DA_PER_BLOCK_UM="${TARGET_DA_PER_BLOCK_UM:-5}"
DA_PHYS="${DA_PHYS:-5e-6}"
TARGET_DB="${TARGET_DB:-0.01}"
TARGET_DN_STORE="${TARGET_DN_STORE:-0.01}"
TARGET_DN_EMIT="${TARGET_DN_EMIT:-0.20}"
TARGET_DN_MOBILE="${TARGET_DN_MOBILE:-0.20}"
CYCLE_PHASES="${CYCLE_PHASES:-8}"
CASE_FILTER="${CASE_FILTER:-}"
KLIST_OVERRIDE="${KLIST_OVERRIDE:-}"
MAKE_2D_PLOTS="${MAKE_2D_PLOTS:-1}"
PRODUCTION_LEVEL="${PRODUCTION_LEVEL:-coarse}"
CRACK_BACKEND="${CRACK_BACKEND:-sharp_wake}"
CZM_MAX_ANGLE_ERROR_DEG="${CZM_MAX_ANGLE_ERROR_DEG:-35}"

should_run_case() {
  local label="$1"
  if [[ -z "${CASE_FILTER}" ]]; then
    return 0
  fi
  local normalized=" ${CASE_FILTER//,/ } "
  [[ "${normalized}" == *" ${label} "* ]]
}

klist_for_level() {
  local coarse="$1"
  local full="$2"
  if [[ -n "${KLIST_OVERRIDE}" ]]; then
    echo "${KLIST_OVERRIDE}"
  elif [[ "${PRODUCTION_LEVEL}" == "full" ]]; then
    echo "${full}"
  else
    echo "${coarse}"
  fi
}

run_case() {
  local label="$1"
  local coarse_klist="$2"
  local full_klist="$3"
  local G00="$4"
  local sigc="$5"
  local expa="$6"
  local expn="$7"
  local floor="$8"
  local emitS="$9"
  local peierlsS="${10}"
  local taylorS="${11}"

  if ! should_run_case "${label}"; then
    echo "=== Skipping ${label} due to CASE_FILTER ==="
    return 0
  fi

  local klist
  klist="$(klist_for_level "${coarse_klist}" "${full_klist}")"
  local plot_flag="--make-2d-plots"
  if [[ "${MAKE_2D_PLOTS}" == "0" || "${MAKE_2D_PLOTS}" == "false" ]]; then
    plot_flag="--no-make-2d-plots"
  fi

  echo "=== Production 2-D run: ${label}; Kmax=[${klist}] ==="
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
    --crack-backend "${CRACK_BACKEND}" \
    --czm-max-angle-error-deg "${CZM_MAX_ANGLE_ERROR_DEG}" \
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

# Coarse lists cover low/no-growth, Paris growth, and high-K rapid fracture.
# Full lists add spacing around cutoffs and Paris slopes. R=0.1, so DeltaK=0.9*Kmax.
run_case "FCC_like_case29" \
  "3.5 4.5 5.5 6.5 7.5 8.5 10.0 12.0" \
  "3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 10.0 11.0 12.0" \
  1.0 2.5 0.70 0.6 0.020 0.375 0.001875 0.0075

run_case "shifted_ductile_case64" \
  "4.0 5.0 6.0 7.0 8.0 9.5 11.0 12.5" \
  "3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 10.0 11.0 12.0 12.5" \
  1.0 3.0 0.70 0.6 0.010 0.375 0.001875 0.0075

run_case "steep_cleavage_case35" \
  "4.0 5.0 6.0 7.0 8.0 9.0 10.0" \
  "3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0" \
  1.0 2.5 0.70 1.0 0.020 0.0 0.0 0.0

run_case "slow_threshold_case101" \
  "5.0 6.0 7.0 8.0 9.5 11.0 12.5" \
  "4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 10.0 11.0 12.0 12.5" \
  1.0 3.5 0.70 0.6 0.020 0.0 0.0 0.0

run_case "higher_barrier_case171" \
  "5.0 6.0 7.0 8.0 9.5 11.0 12.5" \
  "4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 10.0 11.0 12.0 12.5" \
  1.1 2.5 0.70 0.6 0.005 0.0 0.0 0.0

run_case "plastic_shielded_case64_M1" \
  "4.0 5.0 6.0 6.5 7.0 7.5 8.5 10.0 12.0" \
  "3.5 4.0 4.5 5.0 5.5 6.0 6.25 6.5 6.75 7.0 7.25 7.5 8.0 8.5 9.0 10.0 11.0 12.0" \
  1.0 3.0 0.70 0.6 0.010 0.75 0.00375 0.015

python postprocess_v8_material_response_atlas.py \
  --root "${OUTROOT}" \
  --case-table selected_2d_material_response_cases.csv \
  --R "${R}" \
  --cycles-max "${CYCLES_MAX}" \
  --target-crack-extension-um "${TARGET_EXT_UM}" \
  --x-deltaK initial \
  --extract-local-points
