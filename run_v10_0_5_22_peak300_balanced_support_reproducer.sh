#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
FAMILY_JSON=${FAMILY_JSON:-$ROOT/runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json}
OUT=${OUT:-$ROOT/runs/v10_0_5_22_peak300_balanced_support_100um_v1}
TARGET_EXT_UM=${TARGET_EXT_UM:-100}
STEPS=${STEPS:-100000}

OPTION=v913_paper_peak01_0242980_persistent_sites
TEMPERATURE=300
HAZARD_SEED=2085178550

if [[ ! -f "$FAMILY_JSON" ]]; then
  echo "ERROR: missing signed-kernel family: $FAMILY_JSON" >&2
  exit 1
fi

rm -rf "$OUT"
mkdir -p "$OUT"

env \
  CLEAVAGE_HAZARD_MODE=exponential \
  CLEAVAGE_HAZARD_SEED="$HAZARD_SEED" \
  CLEAVAGE_HAZARD_MIN_THRESHOLD=1e-12 \
  CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled \
  CLEAVAGE_EVENT_MIN_FACTOR=0.5 \
  CLEAVAGE_EVENT_MAX_FACTOR=4.0 \
  ARRHENIUS_QUALITY_SUBDIVISION_COUNTS=2,3,4,6,8 \
  ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH=9 \
  ARRHENIUS_QUALITY_PARTITION_MIN_SEGMENT_M=2e-8 \
  python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_22_four_class_standalone \
    --parameter-entry v10.0.5.17-frozen-four-class \
    --parameter-option "$OPTION" \
    --signed-kernel-family "$FAMILY_JSON" \
    --tip-refinement-radius-um 330 \
    --selected-cluster-J-outer-um 240 \
    --local-J-outer-um 100 \
    --mode 2d \
    --bulk-plasticity-mode tip_only \
    --temperatures "$TEMPERATURE" \
    --steps "$STEPS" \
    --nx 36 --ny 72 \
    --tip-h-fine 2.5e-6 --tip-ratio 1.15 \
    --dU 2e-5 --dt 840 \
    --n-stagger 1 \
    --print-every 50 \
    --adaptive-events \
    --adaptive-event-target 0.05 \
    --adaptive-min-frac 1e-8 \
    --adaptive-grow 4 \
    --da-phys 5e-6 \
    --target-crack-extension-um "$TARGET_EXT_UM" \
    --crystal-aniso \
    --crystal-compete \
    --crystal-theta-deg 30 \
    --crystal-C11 523e9 \
    --crystal-C12 203e9 \
    --crystal-C44 160e9 \
    --cleave-gamma-aniso 0.3 \
    --crystal-material w \
    --plane-gate-global \
    --min-global-forward 0.05 \
    --max-fronts 1 \
    --crack-backend adaptive_czm \
    --czm-max-angle-error-deg 35 \
    --j-decomposition cluster \
    --mpz-length-um 50 \
    --mpz-n-bins 80 \
    --save-snapshots 0 \
    --no-plots \
    --out "$OUT" \
  2>&1 | tee "$OUT/console.log"
