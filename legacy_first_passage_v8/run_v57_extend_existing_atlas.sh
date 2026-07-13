#!/usr/bin/env bash
set -euo pipefail

# V5.7 corrected extension runner.
# Run from the Fatigue-PF project root after copying:
#   - this .sh and run_v57_extend_existing_atlas.py into pwd
#   - run_v1_two_barrier_dbtt_fatigue_map_corrected.py into pwd
#   - run_adaptive_two_barrier_threshold_study.py into pwd
#   - fracture_contexts_v5_7.csv into pwd
# with the arrhenius_fracture package importable from pwd.
#
# Stage 1 (vectorized): monotonic Kc(T) per (surface, context, T), resumable.
# Stage 2 (serial V1):  da/dN(DeltaK,T) points and rate-defined DeltaK_th at
#                       1e-10 (primary) and 1e-12 (sensitivity) m/cycle,
#                       resumable at the raw (surface, context, T, DeltaK) level.
#
# Environment toggles:
#   RESUME=1          resume both stages from existing raw CSVs
#   ANALYSIS_ONLY=1   rebuild descriptors/thresholds/associations only
#   SKIP_FATIGUE=1    stage 1 + analysis only
#   FATIGUE_SURFACES="EXP_00012 EXP_00345"   explicit surface list for stage 2

for f in \
  run_v57_extend_existing_atlas.py \
  run_v1_two_barrier_dbtt_fatigue_map_corrected.py \
  run_adaptive_two_barrier_threshold_study.py \
  "${CONTEXT_TABLE:-fracture_contexts_v5_7.csv}" \
  arrhenius_fracture/sharp_front.py \
  arrhenius_fracture/fatigue_v1.py \
  arrhenius_fracture/sn_arrhenius_chain.py \
  arrhenius_fracture/config.py; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing required V5.7 file: $f" >&2
    exit 2
  fi
done

echo "Environment: ${CONDA_DEFAULT_ENV:-none}"
echo "Python: $(command -v python)"

python - <<'PY'
from arrhenius_fracture.config import ElasticProperties, FractureBarrier
from arrhenius_fracture.sharp_front import FrontConfig, FrontEngine
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.fatigue_v1 import (FatigueWaveform, FatigueControllerConfig,
                                           FatigueCycleHazardController)
import importlib.util
for path, req in [
    ("run_v1_two_barrier_dbtt_fatigue_map_corrected.py",
     ["_map_cycle_step", "AnchoredCleavageBarrier", "_compat_barrier_diagnostics"]),
    ("run_adaptive_two_barrier_threshold_study.py",
     ["effective_rate", "locate_bracket", "crossing_estimate", "threshold_record"]),
]:
    spec = importlib.util.spec_from_file_location("preflight_" + path.replace(".", "_"), path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    missing = [r for r in req if not hasattr(mod, r)]
    assert not missing, f"{path} missing {missing}"
print("V5.7 corrected import preflight OK")
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
SKIP_FATIGUE_ARG=""
if [[ "${SKIP_FATIGUE:-0}" == "1" || "${SKIP_FATIGUE:-0}" == "true" ]]; then
  SKIP_FATIGUE_ARG="--skip-fatigue-growth"
fi
FSURF_ARG=""
if [[ -n "${FATIGUE_SURFACES:-}" ]]; then
  FSURF_ARG="--fatigue-surface-ids ${FATIGUE_SURFACES}"
fi
CONTEXT_FILTER_ARG=""
if [[ -n "${CONTEXT_FILTER:-}" ]]; then
  CONTEXT_FILTER_ARG="--context-filter ${CONTEXT_FILTER}"
fi

python run_v57_extend_existing_atlas.py \
  --atlas-dir "${ATLAS_DIR:-runs/sn_v1_barrier_phenomena_map_v5_6}" \
  --out "${OUT:-runs/v5_7_extension}" \
  --context-table "${CONTEXT_TABLE:-fracture_contexts_v5_7.csv}" \
  --temperatures ${TEMPS:-100 200 300 400 500 600 700 800 900} \
  --surface-start "${SURFACE_START:-0}" \
  ${SURFACE_STOP:+--surface-stop ${SURFACE_STOP}} \
  --monotonic-Kmax-MPa "${MONOTONIC_KMAX:-40.0}" \
  --monotonic-dK-MPa "${MONOTONIC_DK:-0.10}" \
  --Kdot-MPa-sqrtm-per-s "${KDOT:-0.005}" \
  --fatigue-temperatures ${FATIGUE_TEMPS:-100 300 500 700 900} \
  --n-fatigue-surfaces "${N_FATIGUE_SURFACES:-96}" \
  --fatigue-surface-seed "${FATIGUE_SURFACE_SEED:-42}" \
  --rate-criteria ${RATE_CRITERIA:-1e-10 1e-12} \
  --primary-rate-criterion "${PRIMARY_RATE_CRITERION:-1e-10}" \
  --DeltaK-seeds ${DELTAK_SEEDS:-0.05 0.10 0.20 0.40 0.80 1.60 3.20 6.40 12.80} \
  --DeltaK-min "${DELTAK_MIN:-0.025}" \
  --DeltaK-max "${DELTAK_MAX:-20.0}" \
  --threshold-abs-tol "${THRESHOLD_ABS_TOL:-0.05}" \
  --threshold-rel-tol "${THRESHOLD_REL_TOL:-0.03}" \
  --max-refine-iters "${MAX_REFINE_ITERS:-10}" \
  --R "${R:-0.1}" \
  --frequency-Hz "${FREQUENCY_HZ:-1000}" \
  --cycles-max "${CYCLES_MAX:-2e14}" \
  --max-blocks "${MAX_BLOCKS:-10000}" \
  --n-advances "${N_ADVANCES:-5}" \
  --da-m "${DA_M:-20e-6}" \
  --n-phase "${N_PHASE:-96}" \
  --target-dB "${TARGET_DB:-0.02}" \
  --target-dN-store "${TARGET_DN_STORE:-0.01}" \
  ${CONTEXT_FILTER_ARG:+$CONTEXT_FILTER_ARG} \
  ${FSURF_ARG:+$FSURF_ARG} \
  ${SKIP_FATIGUE_ARG:+$SKIP_FATIGUE_ARG} \
  ${RESUME_ARG:+$RESUME_ARG} \
  ${ANALYSIS_ARG:+$ANALYSIS_ARG}
