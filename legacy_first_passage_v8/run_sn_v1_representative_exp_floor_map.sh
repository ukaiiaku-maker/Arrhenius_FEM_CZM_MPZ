#!/usr/bin/env bash
set -euo pipefail

# V5.6 manual-overlay runner.
# Run from the Fatigue-PF project root after copying:
#   - this .sh and the matching .py into pwd
#   - the four module files into ./arrhenius_fracture/

for f in \
  run_sn_v1_representative_exp_floor_map.py \
  arrhenius_fracture/sn_arrhenius_chain.py \
  arrhenius_fracture/sn_v1.py \
  arrhenius_fracture/sn_v1_arrhenius.py \
  arrhenius_fracture/sn_v1_arrhenius_batch.py; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing required V5.6 file: $f" >&2
    exit 2
  fi
done

echo "Environment: ${CONDA_DEFAULT_ENV:-none}"
echo "Python: $(command -v python)"

python - <<'PY'
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_v1_arrhenius import SNCase
from arrhenius_fracture.sn_v1_arrhenius_batch import run_stress_grid
print("V5.6 import preflight OK")
PY

python -m compileall -q arrhenius_fracture

RESUME_ARG=""
if [[ "${RESUME:-0}" == "1" || "${RESUME:-0}" == "true" ]]; then
  RESUME_ARG="--resume"
fi
ANALYSIS_ARG=""
if [[ "${ANALYSIS_ONLY:-0}" == "1" || "${ANALYSIS_ONLY:-0}" == "true" ]]; then
  ANALYSIS_ARG="--analysis-only"
fi
RF_ARG=""
if [[ "${RUN_RF:-0}" == "1" || "${RUN_RF:-0}" == "true" ]]; then
  RF_ARG="--run-random-forest"
fi
DESIGN_ARG=""
if [[ -n "${DESIGN_CSV:-}" ]]; then
  DESIGN_ARG="--design-csv ${DESIGN_CSV}"
fi

python run_sn_v1_representative_exp_floor_map.py \
  --out "${OUT:-runs/sn_v1_barrier_phenomena_map_v5_6}" \
  --n-surfaces "${N_SURFACES:-3840}" \
  --candidate-pool "${CANDIDATE_POOL:-2048}" \
  --design-batch-size "${DESIGN_BATCH_SIZE:-256}" \
  --seed "${SEED:-42}" \
  --temperatures ${TEMPS:-25 40 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 450 500 550 600 650 700 750 800 850 900 1000 1100 1200} \
  --prescreen-temperatures ${PRESCREEN_TEMPS:-25 50 75 100 150 200 250 300 400 500 600 750 900 1100 1200} \
  --strength-rate "${STRENGTH_RATE:-1e-4}" \
  --anomaly-gain-min "${ANOMALY_GAIN_MIN:-0.05}" \
  --anomaly-gain-main-max "${ANOMALY_GAIN_MAIN_MAX:-0.50}" \
  --anomaly-gain-sensitivity-max "${ANOMALY_GAIN_SENSITIVITY_MAX:-1.00}" \
  --association-anomaly-caps ${ASSOCIATION_ANOMALY_CAPS:-0.25 0.50 1.00} \
  --fatigue-temperatures ${FATIGUE_TEMPS:-100 200 300 400 500 600 700} \
  --stress-reference-temperature "${STRESS_REFERENCE_T:-300}" \
  --cycles-max "${CYCLES_MAX:-1e12}" \
  --max-blocks "${MAX_BLOCKS:-12000}" \
  --block-cycles "${BLOCK_CYCLES:-1e10}" \
  --n-phase "${N_PHASE:-32}" \
  --target-dP "${TARGET_DP:-0.03}" \
  --target-dD "${TARGET_DD:-0.03}" \
  --target-rho-rel-block "${TARGET_RHO_REL_BLOCK:-0.15}" \
  --target-dB-nuc "${TARGET_DB_NUC:-0.20}" \
  --stress-fractions ${STRESS_FRACTIONS:-0.025 0.05 0.10 0.18 0.30 0.48 0.72 1.00 1.35} \
  --refine-target-lives ${REFINE_TARGET_LIVES:-1e6 1e8 1e10 1e12} \
  --refine-rounds "${REFINE_ROUNDS:-1}" \
  --min-stress-MPa "${MIN_STRESS_MPA:-5}" \
  --max-stress-MPa "${MAX_STRESS_MPA:-5000}" \
  --checkpoint-every "${CHECKPOINT_EVERY:-10}" \
  ${DESIGN_ARG:+$DESIGN_ARG} \
  ${RESUME_ARG:+$RESUME_ARG} \
  ${ANALYSIS_ARG:+$ANALYSIS_ARG} \
  ${RF_ARG:+$RF_ARG}
