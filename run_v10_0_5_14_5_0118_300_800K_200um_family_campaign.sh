#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PFROOT=${PFROOT:-/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1}
FAMILY_JSON=${FAMILY_JSON:-$PFROOT/runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_14_5_0118_300_800K_200um_family_v1}
TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800"}
TARGET_EXT_UM=${TARGET_EXT_UM:-200}
STEPS=${STEPS:-100000}
DU=${DU:-2e-5}
DT=${DT:-840}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-10}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-20}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
PRINT_EVERY=${PRINT_EVERY:-100}

if [[ ! -f "$FAMILY_JSON" ]]; then
  echo "ERROR: missing PF signed-kernel family: $FAMILY_JSON" >&2
  exit 1
fi
mkdir -p "$CAMPAIGN_ROOT"

cat > "$CAMPAIGN_ROOT/campaign_configuration.txt" <<EOF
release=10.0.5.14.5
candidate=v912_peak_0118_persistent_sites
temperatures_K=$TEMPERATURES
target_extension_um=$TARGET_EXT_UM
steps=$STEPS
dU=$DU
dt_s=$DT
max_jobs=$MAX_JOBS
snapshot_spacing_um=$SNAPSHOT_BY_EXT_UM
kernel_family=$FAMILY_JSON
transport=hybrid_asymptotic_retained_exact__physical_BE
adaptive_czm_tip_support=authoritative_plus_minus_pair
EOF

run_case() {
  set -euo pipefail
  local T="$1"
  local TAG
  printf -v TAG "%04d" "$T"
  local OUT="$CAMPAIGN_ROOT/T${TAG}K"
  local MANIFEST="$OUT/persistent_site_production_manifest_v10_0_5_14_5.json"
  local LOG="$OUT/console.log"

  if [[ "$SKIP_FINISHED" == "1" && -f "$MANIFEST" ]] &&
     grep -q '"run_completed_without_exception": true' "$MANIFEST"; then
    echo "[SKIP] T=${T} K already complete"
    return 0
  fi

  rm -rf "$OUT"
  mkdir -p "$OUT"
  echo "============================================================"
  echo "[START] T=${T} K"
  echo "  output: $OUT"
  echo "  target: ${TARGET_EXT_UM} um"
  echo "============================================================"

  python -m \
    arrhenius_fracture.mode_i_first_passage_v10_0_5_14_5_persistent_site_family \
    --persistent-site-option v912_peak_0118_persistent_sites \
    --signed-kernel-family "$FAMILY_JSON" \
    --tip-refinement-radius-um 330 \
    --selected-cluster-J-outer-um 240 \
    --local-J-outer-um 100 \
    --mode 2d \
    --bulk-plasticity-mode tip_only \
    --temperatures "$T" \
    --steps "$STEPS" \
    --nx 36 --ny 72 \
    --tip-h-fine 2.5e-6 --tip-ratio 1.15 \
    --dU "$DU" --dt "$DT" \
    --n-stagger 1 \
    --print-every "$PRINT_EVERY" \
    --adaptive-events \
    --adaptive-event-target 0.05 \
    --adaptive-min-frac 1e-8 \
    --adaptive-grow 4 \
    --da-phys 5e-6 \
    --target-crack-extension-um "$TARGET_EXT_UM" \
    --crystal-aniso \
    --crystal-compete \
    --crystal-theta-deg 45 \
    --crystal-C11 523e9 \
    --crystal-C12 203e9 \
    --crystal-C44 160e9 \
    --cleave-gamma-aniso 0.3 \
    --crystal-material w \
    --max-fronts 1 \
    --crack-backend adaptive_czm \
    --czm-max-angle-error-deg 35 \
    --j-decomposition cluster \
    --mpz-length-um 50 \
    --mpz-n-bins 80 \
    --save-snapshots "$SAVE_SNAPSHOTS" \
    --snapshot-cols "$SNAPSHOT_COLS" \
    --snapshot-by-crack-extension-um "$SNAPSHOT_BY_EXT_UM" \
    --no-plots \
    --out "$OUT" \
    2>&1 | tee "$LOG"
  echo "[DONE] T=${T} K"
}

export -f run_case
export ROOT PFROOT FAMILY_JSON CAMPAIGN_ROOT TARGET_EXT_UM STEPS DU DT
export SKIP_FINISHED SNAPSHOT_BY_EXT_UM SAVE_SNAPSHOTS SNAPSHOT_COLS PRINT_EVERY

printf '%s\n' $TEMPERATURES |
  xargs -n 1 -P "$MAX_JOBS" bash -c 'run_case "$1"' _

echo "Campaign finished: $CAMPAIGN_ROOT"
