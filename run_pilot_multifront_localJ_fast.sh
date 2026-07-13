#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-pilot_multifront_localJ_fast}
mkdir -p "$OUT"
COMMON="--mode 2d --nx 40 --ny 80 \
  --tip-h-fine 0.6e-6 --tip-ratio 1.25 \
  --dU 2e-7 --dt 8.4 --steps 700 --n-stagger 2 \
  --save-snapshots 6 --snapshot-cols 3 --print-every 100 \
  --crystal-aniso --crystal-compete --crystal-branch --crystal-material branchy \
  --cleave-gamma-aniso 2.0 --branch-overdrive-ratio 0.70 \
  --branch-share-mode hazard --branch-hazard-sharpness 2.0 \
  --branch-energy-share hazard-budget \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --emb-sat-frac 1 --n-sat 2000 \
  --adaptive-events --adaptive-event-target 0.35 --adaptive-min-frac 1e-8 \
  --adaptive-grow 4.0 --max-fronts 24 --branch-spacing 10 \
  --rJ-outer 6e-6 --branch-resolve-length 6e-6"
for T in 700 900; do
  for TH in 30 45; do
    CASE="$OUT/T${T}_th${TH}"
    echo "=== $CASE ==="
    python3 -m arrhenius_fracture.sharp_front $COMMON \
      --temperatures "$T" --crystal-theta-deg "$TH" --out "$CASE"
  done
done
