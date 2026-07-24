#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PFROOT=${PFROOT:-/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1}
FAMILY_JSON=${FAMILY_JSON:-$PFROOT/runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json}
OUTROOT=${OUTROOT:-$ROOT/runs/v10_0_5_14_2_0118_700K_20um_family_smoke_v1}

T_K=${T_K:-700}
TARGET_EXT_UM=${TARGET_EXT_UM:-20}
STEPS=${STEPS:-20000}
DU=${DU:-2e-5}
DT=${DT:-840}
PRINT_EVERY=${PRINT_EVERY:-100}

if [[ ! -f "$FAMILY_JSON" ]]; then
  echo "ERROR: missing PF kernel family: $FAMILY_JSON" >&2
  exit 1
fi

mkdir -p "$OUTROOT"

cat <<EOF
v10.0.5.14.2 candidate-0118 kernel-family smoke
  root:        $ROOT
  family:      $FAMILY_JSON
  output:      $OUTROOT
  temperature: ${T_K} K
  target:      ${TARGET_EXT_UM} um
  dU / dt:     $DU / $DT s
  transport:   adaptive implicit backward-Euler upwind
EOF

set -o pipefail
python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_14_2_persistent_site_family \
  --persistent-site-option v912_peak_0118_persistent_sites \
  --signed-kernel-family "$FAMILY_JSON" \
  --tip-refinement-radius-um 330 \
  --selected-cluster-J-outer-um 240 \
  --local-J-outer-um 100 \
  --mode 2d \
  --bulk-plasticity-mode tip_only \
  --temperatures "$T_K" \
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
  --save-snapshots 2 \
  --snapshot-cols 2 \
  --snapshot-by-crack-extension-um 10 \
  --no-plots \
  --out "$OUTROOT" \
  2>&1 | tee "$OUTROOT/console.log"
