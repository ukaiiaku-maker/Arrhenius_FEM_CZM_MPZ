#!/usr/bin/env bash
set -euo pipefail

# 2-D material-response atlas for the v8 sharp-front fatigue model.
# See README_2D_MATERIAL_RESPONSE_ATLAS.md for the case definitions and interpretation.

python -m compileall -q arrhenius_fracture

OUTROOT="${OUTROOT:-runs/v8_material_response_atlas_2d}"
KLIST="${KLIST:-4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0}"
T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
CYCLES_MAX="${CYCLES_MAX:-1e11}"
BLOCKS="${BLOCKS:-220}"
NX="${NX:-24}"
NY="${NY:-48}"
SNAPSHOTS="${SNAPSHOTS:-12}"
SNAPSHOT_COLS="${SNAPSHOT_COLS:-4}"

COMMON="\
  --Kmax-MPa-sqrt-m ${KLIST} \
  --T ${T} \
  --R ${R} \
  --frequency-Hz ${FREQ} \
  --blocks ${BLOCKS} \
  --cycles-max ${CYCLES_MAX} \
  --block-cycles 1e5 \
  --max-block-cycles inf \
  --cycle-block-mode hazard_limited \
  --target-dB 0.02 \
  --target-dN-store 0.025 \
  --target-dN-emit 0.25 \
  --target-dN-mobile 0.25 \
  --storage-model escape_limited \
  --calibrate-2d-K \
  --K-calib-iters 3 \
  --K-calib-tol 5e-3 \
  --no-stop-after-first-2d-fire \
  --cyclic-mechanics-phases 8 \
  --nx ${NX} --ny ${NY} \
  --tip-h-fine 1e-6 \
  --tip-ratio 1.25 \
  --save-snapshots ${SNAPSHOTS} \
  --snapshot-cols ${SNAPSHOT_COLS} \
  --make-2d-plots \
  --min-global-forward 0.05 \
  --cleave-barrier-kind exp_floor \
  --cleave-exp-T-mode mu_scale \
  --emit-energy-scale 0.75 \
  --peierls-energy-scale 0.00375 \
  --taylor-energy-scale 0.015 \
  --peierls-stress-scale 1.0 \
  --taylor-stress-scale 1.0"

run_case() {
  local label="$1"
  local G00="$2"
  local sigc="$3"
  local expa="$4"
  local expn="$5"
  local floor="$6"
  local emitS="$7"
  local peierlsS="$8"
  local taylorS="$9"

  echo "=== Running ${label} ==="
  python run_v8_compare_1d_2d_K_sweep.py ${COMMON} \
    --out "${OUTROOT}/${label}" \
    --cleave-G00-eV "${G00}" \
    --cleave-sigc0-GPa "${sigc}" \
    --cleave-exp-a "${expa}" \
    --cleave-exp-n "${expn}" \
    --cleave-floor-frac "${floor}" \
    --emit-entropy-scale "${emitS}" \
    --peierls-entropy-scale "${peierlsS}" \
    --taylor-entropy-scale "${taylorS}"
}

# E = 0.75, rP = 0.005, rT = 0.02.
# M = 0.5 -> emitS=0.375, peierlsS=0.001875, taylorS=0.0075.
# M = 0.0 -> entropy scales all zero.
# M = 1.0 -> emitS=0.75, peierlsS=0.00375, taylorS=0.015.

run_case "FCC_like_case29"             1.0 2.5 0.70 0.6 0.020 0.375 0.001875 0.0075
run_case "shifted_ductile_case64"      1.0 3.0 0.70 0.6 0.010 0.375 0.001875 0.0075
run_case "steep_cleavage_case35"       1.0 2.5 0.70 1.0 0.020 0.0   0.0      0.0
run_case "slow_threshold_case101"      1.0 3.5 0.70 0.6 0.020 0.0   0.0      0.0
run_case "higher_barrier_case171"      1.1 2.5 0.70 0.6 0.005 0.0   0.0      0.0
run_case "plastic_shielded_case64_M1"  1.0 3.0 0.70 0.6 0.010 0.75  0.00375  0.015

python postprocess_v8_material_response_atlas.py \
  --root "${OUTROOT}" \
  --case-table selected_2d_material_response_cases.csv \
  --R "${R}"
