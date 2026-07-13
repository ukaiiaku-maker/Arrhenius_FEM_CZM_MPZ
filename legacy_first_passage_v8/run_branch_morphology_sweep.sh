#!/usr/bin/env bash
set -euo pipefail

# Morphology-resolved sweep for branch path shapes.
# Purpose: avoid the one-step runaway seen in the smoke test by using small load/time
# increments. Do NOT use this for rate-calibrated Kc/DBTT; use the birth sweep for that.

ROOT=${1:-branch_morphology_sweep_$(date +%Y%m%d_%H%M)}
mkdir -p "$ROOT"

COMMON=(
  --mode 2d
  --nx 50 --ny 100
  --tip-h-fine 0.6e-6 --tip-ratio 1.25
  --n-stagger 2
  --save-snapshots 8 --snapshot-cols 4
  --print-every 25
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
  --v-rayleigh 2600
)

# Smaller axis set: these are the most branch-informative orientations.
TEMPS=(700 850 900 950 1100)
THETAS=(30 37.5 45 52.5 60)
GAMMAS=(1.0 2.0 3.0)
RATIOS=(0.55 0.70 0.85)

# Small increments to observe path evolution. This is morphology-only, not rate-matched.
DU=1e-7
DT=1e-7
STEPS=260

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

python3 analyze_branch_sweep.py "$ROOT" --out "$ROOT/branch_sweep_summary.csv"
echo "Done. Summary: $ROOT/branch_sweep_summary.csv"
