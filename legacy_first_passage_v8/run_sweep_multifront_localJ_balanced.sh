#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-sweep_multifront_localJ_balanced}
mkdir -p "$OUT"

# Balanced local-J recipe:
#   da_phys < rJ_outer < branch_resolve_length
#   branch_resolve_length is several increments, so newborn branches co-grow
#   under the parent budget until they are longer than the J/process-zone scale.
COMMON="--mode 2d --nx 50 --ny 100 \
  --tip-h-fine 0.30e-6 --tip-ratio 1.20 \
  --dU 2e-7 --dt 8.4 --steps 2500 --n-stagger 2 \
  --save-snapshots 8 --snapshot-cols 4 --print-every 200 \
  --crystal-aniso --crystal-compete --crystal-branch --crystal-material branchy \
  --cleave-gamma-aniso 2.0 --branch-overdrive-ratio 0.70 \
  --branch-share-mode hazard --branch-hazard-sharpness 2.0 \
  --branch-energy-share hazard-budget \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --emb-sat-frac 1 --n-sat 2000 \
  --adaptive-events --adaptive-event-target 0.35 --adaptive-min-frac 1e-8 \
  --adaptive-grow 4.0 --max-fronts 32 \
  --da-phys 2e-6 --branch-spacing 15 \
  --rJ-outer 12e-6 --branch-resolve-length 20e-6"

for T in 700 900 1100; do
  for TH in 30 45; do
    CASE="$OUT/T${T}_th${TH}"
    echo "=== $CASE ==="
    python3 -m arrhenius_fracture.sharp_front $COMMON \
      --temperatures "$T" --crystal-theta-deg "$TH" --out "$CASE"
  done
done
python3 "$(dirname "$0")/analyze_localJ_growth.py" "$OUT" --out "$OUT/localJ_growth_summary.csv" || true
