#!/usr/bin/env bash
set -euo pipefail

# 1-D temperature sweep for the six material-response classes used in the
# production 2-D atlas.  This uses the same V1 fatigue_sharp_front engine and
# the same cleavage/plasticity settings as the production comparison runs.
#
# Main modes:
#   ENTROPY_MODE=as_calibrated  : preserve the 300 K atlas settings exactly
#   ENTROPY_MODE=common_m30     : set emission S* ~ -30 kB in all six cases
#   ENTROPY_MODE=common_m40     : set emission S* ~ -40 kB in all six cases
#   ENTROPY_MODE=common_m50     : set emission S* ~ -50 kB in all six cases
#
# For the common physical-entropy modes, Peierls/Taylor entropy scales are set
# to zero so the temperature sensitivity isolates crack-tip emission versus
# crack opening. Their energy barriers remain unchanged.
#
# Useful overrides:
#   CASE_FILTER="plastic_shielded_case64_M1"
#   TEMPS="300 400 500 600 700 900"
#   KLIST="3.5 4.5 5.5 6.5 7.5 8.5 10 12 14"
#   CYCLES_MAX=2e14 MAX_BLOCKS=5000
#   ENTROPY_MODE=common_m40

python -m compileall -q arrhenius_fracture

OUTROOT="${OUTROOT:-runs/v1_material_response_temperature_sweep}"
TEMPS="${TEMPS:-300 400 500 600 700 900}"
KLIST="${KLIST:-3.5 4.5 5.5 6.5 7.5 8.5 10 12 14}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
CYCLES_MAX="${CYCLES_MAX:-2e14}"
MAX_BLOCKS="${MAX_BLOCKS:-5000}"
ENTROPY_MODE="${ENTROPY_MODE:-as_calibrated}"
CASE_FILTER="${CASE_FILTER:-}"
N_ADVANCES="${N_ADVANCES:-5}"

# W[100] surrogate: gT=0.003934 eV/K.  Because entropy_scale multiplies gT
# independently, the direct scales corresponding to constant low-stress
# entropies are approximately:
#   -30 kB -> 0.657143
#   -40 kB -> 0.876190
#   -50 kB -> 1.095238
entropy_scale_for_mode() {
  case "$ENTROPY_MODE" in
    common_m30) echo "0.657143" ;;
    common_m40) echo "0.876190" ;;
    common_m50) echo "1.095238" ;;
    as_calibrated) echo "" ;;
    *) echo "ERROR: unknown ENTROPY_MODE=${ENTROPY_MODE}" >&2; exit 2 ;;
  esac
}

should_run_case() {
  local label="$1"
  if [[ -z "$CASE_FILTER" ]]; then
    return 0
  fi
  local normalized=" ${CASE_FILTER//,/ } "
  [[ "$normalized" == *" ${label} "* ]]
}

run_case() {
  local label="$1"
  local G00="$2"
  local sigc="$3"
  local expa="$4"
  local expn="$5"
  local floor="$6"
  local emitS_cal="$7"
  local peierlsS_cal="$8"
  local taylorS_cal="$9"

  if ! should_run_case "$label"; then
    echo "=== Skipping ${label} due to CASE_FILTER ==="
    return 0
  fi

  local emitS="$emitS_cal"
  local peierlsS="$peierlsS_cal"
  local taylorS="$taylorS_cal"
  local commonS
  commonS="$(entropy_scale_for_mode)"
  if [[ -n "$commonS" ]]; then
    emitS="$commonS"
    peierlsS="0.0"
    taylorS="0.0"
  fi

  echo "=== V1 T sweep: ${label}; entropy_mode=${ENTROPY_MODE}; emit_entropy_scale=${emitS} ==="

  for K in $KLIST; do
    local KTAG
    KTAG="${K/./p}"
    local OUT="${OUTROOT}/${ENTROPY_MODE}/${label}/K${KTAG}"

    python -m arrhenius_fracture.fatigue_sharp_front \
      --temperatures ${TEMPS} \
      --Kmax-MPa-sqrt-m "$K" \
      --cycles-max "$CYCLES_MAX" \
      --max-blocks "$MAX_BLOCKS" \
      --no-plots \
      --out "$OUT" \
      --R "$R" \
      --frequency-Hz "$FREQ" \
      --block-cycles 1e5 \
      --min-block-cycles 1e-6 \
      --max-block-cycles inf \
      --cycle-block-mode hazard_limited \
      --target-dB 0.02 \
      --target-dN-store 0.01 \
      --target-dN-emit 0.20 \
      --target-dN-mobile 0.20 \
      --target-dN-escape inf \
      --target-dN-peierls inf \
      --target-dN-taylor inf \
      --storage-model escape_limited \
      --dN-cap inf \
      --continue-after-fire \
      --n-advances "$N_ADVANCES" \
      --cleave-barrier-kind exp_floor \
      --cleave-G00-eV "$G00" \
      --cleave-sigc0-GPa "$sigc" \
      --cleave-exp-a "$expa" \
      --cleave-exp-n "$expn" \
      --cleave-floor-frac "$floor" \
      --cleave-exp-T-mode mu_scale \
      --emit-energy-scale 0.75 \
      --emit-entropy-scale "$emitS" \
      --peierls-energy-scale 0.00375 \
      --peierls-entropy-scale "$peierlsS" \
      --peierls-stress-scale 1.0 \
      --taylor-energy-scale 0.015 \
      --taylor-entropy-scale "$taylorS" \
      --taylor-stress-scale 1.0
  done
}

# Exact six base cases used in the 2-D atlas.
run_case "FCC_like_case29"             1.0 2.5 0.70 0.6 0.020 0.375 0.001875 0.0075
run_case "shifted_ductile_case64"      1.0 3.0 0.70 0.6 0.010 0.375 0.001875 0.0075
run_case "steep_cleavage_case35"       1.0 2.5 0.70 1.0 0.020 0.0   0.0      0.0
run_case "slow_threshold_case101"      1.0 3.5 0.70 0.6 0.020 0.0   0.0      0.0
run_case "higher_barrier_case171"      1.1 2.5 0.70 0.6 0.005 0.0   0.0      0.0
run_case "plastic_shielded_case64_M1"  1.0 3.0 0.70 0.6 0.010 0.75  0.00375  0.015

python postprocess_v1_material_response_temperature.py \
  --root "${OUTROOT}/${ENTROPY_MODE}" \
  --R "$R" \
  --cycles-max "$CYCLES_MAX"
