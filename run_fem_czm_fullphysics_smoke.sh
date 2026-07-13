#!/usr/bin/env bash
set -euo pipefail

# Full-physics migration smoke test.
# Uses the existing v8 fatigue/plasticity/process-zone/branch machinery and
# changes only the 2-D crack geometry backend to the Arrhenius discrete CZM.

OUT="${OUT:-runs/fem_czm_fullphysics_smoke}"
KMAX="${KMAX:-7.0}"
T="${T:-300}"

python run_v8_compare_1d_2d_K_sweep.py \
  --out "${OUT}" \
  --Kmax-MPa-sqrt-m "${KMAX}" \
  --T "${T}" \
  --R 0.1 \
  --frequency-Hz 1000 \
  --blocks 3 \
  --cycles-max 1e4 \
  --block-cycles 1 \
  --max-block-cycles 1 \
  --cycle-block-mode hazard_limited \
  --target-dB 0.1 \
  --target-dN-store 0.1 \
  --storage-model escape_limited \
  --nx 12 --ny 24 \
  --tip-h-fine 5e-6 \
  --tip-ratio 1.30 \
  --da-phys 5e-6 \
  --crack-backend edge_split_czm \
  --czm-event-damage 1.0 \
  --czm-max-angle-error-deg 35 \
  --save-snapshots 0 \
  --no-make-2d-plots \
  --no-stop-after-first-2d-fire \
  --no-calibrate-2d-K \
  --cleave-barrier-kind exp_floor \
  --cleave-exp-T-mode mu_scale \
  --cleave-G00-eV 1.0 \
  --cleave-sigc0-GPa 3.0 \
  --cleave-exp-a 0.70 \
  --cleave-exp-n 0.6 \
  --cleave-floor-frac 0.01 \
  --emit-energy-scale 0.75 \
  --emit-entropy-scale 0.75 \
  --peierls-energy-scale 0.00375 \
  --peierls-entropy-scale 0.00375 \
  --taylor-energy-scale 0.015 \
  --taylor-entropy-scale 0.015
