#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-branch_reduced_matrix_v1}
mkdir -p "$ROOT"

# Reduced matrix for morphology-resolved sharp-front branching.
# Keep loading rate matched; slow morphology using max-advances-per-step, not tiny dt.
COMMON=(
  --mode 2d
  --steps 900 --nx 50 --ny 100
  --tip-h-fine 0.6e-6 --tip-ratio 1.25
  --dU 2e-7 --dt 8.4
  --n-stagger 2
  --save-snapshots 4 --snapshot-cols 4
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

# Stage A: smallest useful map. 9 cases.
TEMPS=(700 900 1100)
THETAS=(30 45 60)
GAMMA=2.0
RATIO=0.70

for T in "${TEMPS[@]}"; do
  for THETA in "${THETAS[@]}"; do
    OUT="$ROOT/A_T${T}_th${THETA}_g${GAMMA}_r${RATIO}"
    echo "[A] T=$T theta=$THETA gamma=$GAMMA ratio=$RATIO -> $OUT"
    python3 -m arrhenius_fracture.sharp_front "${COMMON[@]}" \
      --temperatures "$T" \
      --crystal-theta-deg "$THETA" \
      --cleave-gamma-aniso "$GAMMA" \
      --branch-overdrive-ratio "$RATIO" \
      --out "$OUT"
  done
done

# Stage B: local sensitivity around the branch-prone orientation. 9 cases.
# This determines whether branching is robust or only a permissive-threshold artifact.
T=900
THETA=45
GAMMAS=(1.5 2.0 2.5)
RATIOS=(0.55 0.70 0.85)

for GAMMA in "${GAMMAS[@]}"; do
  for RATIO in "${RATIOS[@]}"; do
    OUT="$ROOT/B_T${T}_th${THETA}_g${GAMMA}_r${RATIO}"
    echo "[B] T=$T theta=$THETA gamma=$GAMMA ratio=$RATIO -> $OUT"
    python3 -m arrhenius_fracture.sharp_front "${COMMON[@]}" \
      --temperatures "$T" \
      --crystal-theta-deg "$THETA" \
      --cleave-gamma-aniso "$GAMMA" \
      --branch-overdrive-ratio "$RATIO" \
      --out "$OUT"
  done
done

python3 "$(dirname "$0")/summarize_branch_reduced.py" "$ROOT" --out "$ROOT/reduced_branch_summary.csv"
echo "Wrote $ROOT/reduced_branch_summary.csv"
