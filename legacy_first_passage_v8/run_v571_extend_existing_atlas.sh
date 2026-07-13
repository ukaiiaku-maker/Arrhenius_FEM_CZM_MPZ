#!/usr/bin/env bash
set -euo pipefail

ATLAS_DIR="${ATLAS_DIR:-runs/sn_v1_barrier_phenomena_map_v5_6}"
OUT="${OUT:-runs/sn_v1_barrier_phenomena_extension_v5_7_1}"
CONTEXT_TABLE="${CONTEXT_TABLE:-fracture_contexts_v5_7.csv}"
FRACTURE_TEMPS="${FRACTURE_TEMPS:-100 200 300 400 500 600 700 800 900}"
THRESHOLD_TEMPS="${THRESHOLD_TEMPS:-100 200 300 400 500 600 700 800 900}"
SURFACE_START="${SURFACE_START:-0}"
SURFACE_STOP="${SURFACE_STOP:-}"
THRESHOLD_SURFACE_START="${THRESHOLD_SURFACE_START:-${SURFACE_START}}"
THRESHOLD_SURFACE_STOP="${THRESHOLD_SURFACE_STOP:-${SURFACE_STOP}}"
CONTEXT_FILTER="${CONTEXT_FILTER:-}"
THRESHOLD_CONTEXTS="${THRESHOLD_CONTEXTS:-}"
RATE_CRITERIA="${RATE_CRITERIA:-1e-10 1e-12}"
THRESHOLD_CYCLES_MAX="${THRESHOLD_CYCLES_MAX:-2e14}"
THRESHOLD_MAX_BLOCKS="${THRESHOLD_MAX_BLOCKS:-10000}"
THRESHOLD_N_PHASE="${THRESHOLD_N_PHASE:-96}"
THRESHOLD_N_ADVANCES="${THRESHOLD_N_ADVANCES:-5}"
THRESHOLD_DK_MIN="${THRESHOLD_DK_MIN:-0.025}"
THRESHOLD_DK_MAX="${THRESHOLD_DK_MAX:-20}"
THRESHOLD_ABS_TOL="${THRESHOLD_ABS_TOL:-0.05}"
THRESHOLD_REL_TOL="${THRESHOLD_REL_TOL:-0.03}"
SKIP_MONOTONIC="${SKIP_MONOTONIC:-0}"
SKIP_THRESHOLDS="${SKIP_THRESHOLDS:-0}"
ANALYSIS_ONLY="${ANALYSIS_ONLY:-0}"

export PYTHONNOUSERSITE=1

echo "Environment: ${CONDA_DEFAULT_ENV:-none}"
echo "Python: $(command -v python)"
echo "Atlas: ${ATLAS_DIR}"
echo "Extension out: ${OUT}"

for f in \
  run_v571_extend_existing_atlas.py \
  run_v1_two_barrier_dbtt_fatigue_map_corrected.py \
  "${CONTEXT_TABLE}"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing required file: $f" >&2
    exit 2
  fi
done

python - <<'PY'
import run_v571_extend_existing_atlas
import run_v1_two_barrier_dbtt_fatigue_map_corrected
print("V5.7.1 import preflight OK")
PY

ARGS=(
  --atlas-dir "${ATLAS_DIR}"
  --out "${OUT}"
  --context-table "${CONTEXT_TABLE}"
  --temperatures ${FRACTURE_TEMPS}
  --threshold-temperatures ${THRESHOLD_TEMPS}
  --surface-start "${SURFACE_START}"
  --threshold-surface-start "${THRESHOLD_SURFACE_START}"
  --rate-criteria ${RATE_CRITERIA}
  --threshold-cycles-max "${THRESHOLD_CYCLES_MAX}"
  --threshold-max-blocks "${THRESHOLD_MAX_BLOCKS}"
  --threshold-n-phase "${THRESHOLD_N_PHASE}"
  --threshold-n-advances "${THRESHOLD_N_ADVANCES}"
  --threshold-DeltaK-min "${THRESHOLD_DK_MIN}"
  --threshold-DeltaK-max "${THRESHOLD_DK_MAX}"
  --threshold-abs-tol-MPa-sqrtm "${THRESHOLD_ABS_TOL}"
  --threshold-rel-tol "${THRESHOLD_REL_TOL}"
  --resume
)

if [[ -n "${SURFACE_STOP}" ]]; then
  ARGS+=(--surface-stop "${SURFACE_STOP}")
fi
if [[ -n "${THRESHOLD_SURFACE_STOP}" ]]; then
  ARGS+=(--threshold-surface-stop "${THRESHOLD_SURFACE_STOP}")
fi
if [[ -n "${CONTEXT_FILTER}" ]]; then
  ARGS+=(--context-filter ${CONTEXT_FILTER})
fi
if [[ -n "${THRESHOLD_CONTEXTS}" ]]; then
  ARGS+=(--threshold-context-filter ${THRESHOLD_CONTEXTS})
fi
if [[ "${SKIP_MONOTONIC}" == "1" || "${SKIP_MONOTONIC}" == "true" ]]; then
  ARGS+=(--skip-monotonic)
fi
if [[ "${SKIP_THRESHOLDS}" == "1" || "${SKIP_THRESHOLDS}" == "true" ]]; then
  ARGS+=(--skip-thresholds)
fi
if [[ "${ANALYSIS_ONLY}" == "1" || "${ANALYSIS_ONLY}" == "true" ]]; then
  ARGS+=(--analysis-only)
fi

python run_v571_extend_existing_atlas.py "${ARGS[@]}"
