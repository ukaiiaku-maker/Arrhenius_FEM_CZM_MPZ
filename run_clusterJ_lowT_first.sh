#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-clusterJ_lowT_first_v8}
mkdir -p "$OUT"
COMMON="--mode 2d --nx 50 --ny 100 \
  --tip-h-fine 0.35e-6 --tip-ratio 1.25 \
  --dU 2e-7 --dt 8.4 --steps 2500 --n-stagger 2 \
  --save-snapshots 8 --snapshot-cols 4 --print-every 100 \
  --crystal-aniso --crystal-compete --crystal-branch --crystal-material branchy \
  --cleave-gamma-aniso 2.0 --branch-overdrive-ratio 0.95 \
  --branch-fp-min-ratio 0.98 --branch-clock-target 1.0 \
  --branch-clock-angle-tol-deg 10 \
  --branch-secondary-min-K-ratio 0.90 --branch-secondary-min-K-MPa 10.0 \
  --branch-secondary-min-lambda 0.0 \
  --branch-starved-suppression-radius 120e-6 \
  --branch-starved-K-ratio 0.40 --branch-starved-lambda 1e-20 \
  --branch-starved-max-length-factor 4.0 \
  --retire-stagnant-branches --branch-stagnant-lag 120e-6 \
  --branch-stagnant-K-ratio 0.25 --branch-stagnant-lambda 1e-25 \
  --branch-stagnant-steps 25 --branch-stagnant-no-fire-steps 80 \
  --branch-share-mode hazard --branch-hazard-sharpness 2.0 \
  --branch-energy-share hazard-budget \
  --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
  --multihit-m 3 --multihit-tau 1e-6 \
  --cleave-H0-eV 2.6 --cleave-shield-chi 0.2 \
  --emb-sat-frac 1 --n-sat 2000 \
  --adaptive-events --adaptive-event-target 0.35 --adaptive-min-frac 1e-8 \
  --adaptive-grow 4.0 --max-fronts 32 --branch-spacing 50 \
  --da-phys 2e-6 --j-decomposition cluster \
  --rJ-cluster 10e-6 --rJ-outer 12e-6 --branch-resolve-length 20e-6 \
  --local-J-clearance-factor 1.0 --min-J-active-elems 12 \
  --local-J-handoff-min-K-ratio 0.60"
for TH in 30 45; do
  CASE="$OUT/T300_th${TH}"
  echo "=== $CASE ==="
  python3 -m arrhenius_fracture.sharp_front $COMMON \
    --temperatures 300 --crystal-theta-deg "$TH" --out "$CASE"
done
python3 analyze_clusterJ_growth.py "$OUT" > "$OUT/clusterJ_summary.csv" || true
