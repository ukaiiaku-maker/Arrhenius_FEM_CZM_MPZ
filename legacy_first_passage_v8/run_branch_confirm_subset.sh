#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-branch_confirm_subset_v1}
mkdir -p "$ROOT"

# Higher-confidence confirmation for selected cases after the reduced matrix.
# Use more snapshots and a slightly finer mesh; still morphology-capped.
COMMON=(
  --mode 2d
  --steps 1200 --nx 60 --ny 120
  --tip-h-fine 0.5e-6 --tip-ratio 1.18
  --dU 2e-7 --dt 8.4
  --n-stagger 2
  --save-snapshots 8 --snapshot-cols 4
  --print-every 100
  --crystal-aniso --crystal-compete --crystal-branch
  --crystal-material branchy
  --branch-share-mode hazard
  --branch-hazard-sharpness 2.0
  --branch-energy-share hazard-budget
  --emit-S-T-c0-kB=-20
  --emit-S-T-c1=0.02
  --emit-S-sigma-max-kB=8
  --multihit-m 3 --multihit-tau 1e-6
  --cleave-H0-eV 2.6
  --cleave-shield-chi 0.2
  --emb-sat-frac 1
  --n-sat 2000
  --v-rayleigh inf
  --max-advances-per-step 1
  --da-phys 2e-6
)

# Edit this list based on reduced_branch_summary.csv.
# Format: T theta gamma ratio label
CASES=(
  "900 45 2.0 0.70 nominal_branch"
  "900 30 2.0 0.70 biased_branch"
  "900 60 2.0 0.70 alternate_branch"
  "1100 45 2.0 0.70 highT_branch"
)

for CASE in "${CASES[@]}"; do
  read -r T THETA GAMMA RATIO LABEL <<< "$CASE"
  OUT="$ROOT/${LABEL}_T${T}_th${THETA}_g${GAMMA}_r${RATIO}"
  echo "[confirm] $LABEL: T=$T theta=$THETA gamma=$GAMMA ratio=$RATIO -> $OUT"
  python3 -m arrhenius_fracture.sharp_front "${COMMON[@]}" \
    --temperatures "$T" \
    --crystal-theta-deg "$THETA" \
    --cleave-gamma-aniso "$GAMMA" \
    --branch-overdrive-ratio "$RATIO" \
    --out "$OUT"
done

python3 "$(dirname "$0")/summarize_branch_reduced.py" "$ROOT" --out "$ROOT/confirm_branch_summary.csv"
echo "Wrote $ROOT/confirm_branch_summary.csv"
