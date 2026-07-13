#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-branch_morphology_sweep_v2_$(date +%Y%m%d_%H%M)}
mkdir -p "$ROOT"

COMMON=(
  --mode 2d
  --nx 50 --ny 100
  --tip-h-fine 0.6e-6 --tip-ratio 1.25
  --n-stagger 2
  --save-snapshots 8 --snapshot-cols 4
  --print-every 50
  --crystal-aniso --crystal-compete --crystal-branch
  --crystal-material branchy
  --branch-share-mode hazard
  --branch-hazard-sharpness 2.0
  --branch-energy-share hazard-budget
  --emit-S-T-c0-kB=-20
  --emit-S-T-c1=0.02
  --emit-S-sigma-max-kB=8
  --multihit-m 3 --multihit-tau 1e-6
  --emb-sat-frac 1
  --n-sat 2000
  --cleave-shield-chi 0.2
  --cleave-H0-eV 2.6
  --v-rayleigh inf
  --max-advances-per-step 1
  --da-phys 2e-6
)

# Pilot morphology axes. Expand after confirming opening/branching.
TEMPS=(700 850 900 950 1100)
THETAS=(37.5 45 52.5)
GAMMAS=(2.0 3.0)
RATIOS=(0.70 0.85)

# Rate-matched small increments. Previous v1 used DT=1e-7, which killed the
# Arrhenius clock. Here dU/dt matches the K-rate family used in the 2-D Kc sweeps.
DU=2e-7
DT=8.4
STEPS=900

for gamma in "${GAMMAS[@]}"; do
  for ratio in "${RATIOS[@]}"; do
    for theta in "${THETAS[@]}"; do
      for T in "${TEMPS[@]}"; do
        tag="T${T}_g${gamma}_r${ratio}_th${theta}"
        out="$ROOT/$tag"
        echo "=== running $tag ==="
        python3 -m arrhenius_fracture.sharp_front "${COMMON[@]}" \
          --temperatures "$T" \
          --crystal-theta-deg "$theta" \
          --cleave-gamma-aniso "$gamma" \
          --branch-overdrive-ratio "$ratio" \
          --dU "$DU" --dt "$DT" --steps "$STEPS" \
          --out "$out" \
          > "$out.log" 2>&1
      done
    done
  done
done

python3 analyze_branch_sweep_v2.py "$ROOT" --out "$ROOT/branch_sweep_summary.csv"
echo "Done. Summary: $ROOT/branch_sweep_summary.csv"
