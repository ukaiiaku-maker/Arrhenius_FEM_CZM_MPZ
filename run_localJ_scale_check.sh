#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-localJ_scale_check}
mkdir -p "$OUT"

BASE="--mode 2d --nx 50 --ny 100 \
  --tip-h-fine 0.30e-6 --tip-ratio 1.20 \
  --dU 2e-7 --dt 8.4 --steps 1400 --n-stagger 2 \
  --save-snapshots 5 --snapshot-cols 5 --print-every 200 \
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
  --da-phys 2e-6 --branch-spacing 15"

for RJ in 8e-6 12e-6 16e-6; do
  case "$RJ" in
    8e-6)  RES=16e-6 ;;
    12e-6) RES=20e-6 ;;
    16e-6) RES=28e-6 ;;
  esac
  CASE="$OUT/T900_th45_rJ${RJ}_resolve${RES}"
  echo "=== $CASE ==="
  python3 -m arrhenius_fracture.sharp_front $BASE \
    --rJ-outer "$RJ" --branch-resolve-length "$RES" \
    --temperatures 900 --crystal-theta-deg 45 --out "$CASE"
done
python3 "$(dirname "$0")/analyze_localJ_growth.py" "$OUT" --out "$OUT/localJ_scale_summary.csv" || true
