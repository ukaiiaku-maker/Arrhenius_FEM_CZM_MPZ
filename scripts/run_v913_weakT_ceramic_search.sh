#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}}"
PYTHON_BIN="${PYTHON_BIN:-python}"
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"

ANCHOR_REGISTRY="${ANCHOR_REGISTRY:-runs/v9_13_weakT_ceramic_paper_handoff_v1/v9_13_weakT_ceramic_paper_handoff.csv}"
BASE_PHYSICS_JSON="${BASE_PHYSICS_JSON:-mpz_v9_13_v10222_transfer_common_physics.json}"
LOADING_MAP="${LOADING_MAP:-runs/v9_13_long_map_exponential_110um_v2/v10_2_22_long_rcurve_loading_map_exponential_110um.json}"
SEARCH_POLICY="${SEARCH_POLICY:-mpz_v9_13_zero_d_weakT_ceramic_search_policy.json}"
LEGACY_POLICY="${LEGACY_POLICY:-mpz_v9_12_targeted_local_search_policy.json}"
OUTROOT="${OUTROOT:-runs/v9_13_weakT_ceramic_search_v1}"
TEMPERATURES_K="${TEMPERATURES_K:-300 600 900 1100 1200}"
TARGET_EXT_UM="${TARGET_EXT_UM:-100}"
SAMPLES="${SAMPLES:-131072}"
ZERO_D_EXACT_PER_CLASS="${ZERO_D_EXACT_PER_CLASS:-128}"
ZERO_D_PROMOTE_PER_CLASS="${ZERO_D_PROMOTE_PER_CLASS:-8}"
FINAL_PROMOTE_PER_CLASS="${FINAL_PROMOTE_PER_CLASS:-3}"
MAX_JOBS="${MAX_JOBS:-4}"
PROGRESS_INTERVAL_S="${PROGRESS_INTERVAL_S:-60}"

TARGET_TAG="${TARGET_EXT_UM//./p}"
ZERO_D="$OUTROOT/zero_d"
ONE_D="$OUTROOT/one_d_${TARGET_TAG}um"
FINAL_ANALYSIS="$OUTROOT/analysis_${TARGET_TAG}um"

for required in \
  "$ANCHOR_REGISTRY" \
  "$BASE_PHYSICS_JSON" \
  "$LOADING_MAP" \
  "$SEARCH_POLICY" \
  "$LEGACY_POLICY"
do
  test -f "$required" || {
    echo "ERROR: missing required input: $required" >&2
    exit 1
  }
done

mkdir -p "$ZERO_D" "$ONE_D" "$FINAL_ANALYSIS"
read -r -a TEMPERATURE_ARRAY <<< "$TEMPERATURES_K"

if [[ "${#TEMPERATURE_ARRAY[@]}" -lt 5 ]]; then
  echo "ERROR: at least five temperatures are required" >&2
  exit 1
fi

"$PYTHON_BIN" - "$LOADING_MAP" "$TARGET_EXT_UM" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
required = float(sys.argv[2])
payload = json.loads(path.read_text())
coverage = sum(float(value) for value in payload["projected_advances_m"]) * 1.0e6
if coverage + 1.0e-9 < required:
    raise SystemExit(
        f"ERROR: loading map {path} covers {coverage:.6g} um; required {required:g} um"
    )
print(f"LOADING_MAP_OK path={path} coverage_um={coverage:.9g} required_um={required:g}")
PY

"$PYTHON_BIN" -m pytest -q \
  tests/test_v913_weakT_ceramic_objectives.py \
  tests/test_v913_weakT_ceramic_pipeline.py

"$PYTHON_BIN" -u scripts/run_v913_zero_d_weakT_ceramic_search.py \
  --anchor-registry "$ANCHOR_REGISTRY" \
  --base-physics-json "$BASE_PHYSICS_JSON" \
  --loading-map "$LOADING_MAP" \
  --policy-json "$SEARCH_POLICY" \
  --samples "$SAMPLES" \
  --exact-per-class "$ZERO_D_EXACT_PER_CLASS" \
  --promote-per-class "$ZERO_D_PROMOTE_PER_CLASS" \
  --jobs "$MAX_JOBS" \
  --temperatures-K "${TEMPERATURE_ARRAY[@]}" \
  --proxy-extension-um "$TARGET_EXT_UM" \
  --exact-extension-um "$TARGET_EXT_UM" \
  --checkpoint-um "$TARGET_EXT_UM" \
  --progress-interval-s "$PROGRESS_INTERVAL_S" \
  --out "$ZERO_D"

CANDIDATE_COUNT="$((2 * ZERO_D_PROMOTE_PER_CLASS))"
CASE_COUNT="$((CANDIDATE_COUNT * ${#TEMPERATURE_ARRAY[@]}))"
echo "V913_WEAKT_CERAMIC_1D_START candidates=$CANDIDATE_COUNT temperatures=${#TEMPERATURE_ARRAY[@]} cases=$CASE_COUNT target_um=$TARGET_EXT_UM"

SEARCH_STATUS=0
if "$PYTHON_BIN" -u -m scripts.run_v913_autonomous_dbtt_search \
  --candidate-registry "$ZERO_D/combined_promoted_registry.csv" \
  --base-physics-json "$BASE_PHYSICS_JSON" \
  --loading-map "$LOADING_MAP" \
  --policy-json "$LEGACY_POLICY" \
  --families \
  --per-parent 0 \
  --temperatures "${TEMPERATURE_ARRAY[@]}" \
  --checkpoint-um "$TARGET_EXT_UM" \
  --target-extension-um "$TARGET_EXT_UM" \
  --translation-action-exponent 0.95 \
  --max-hazard-increment 0.05 \
  --jobs "$MAX_JOBS" \
  --progress-interval-s "$PROGRESS_INTERVAL_S" \
  --promote-count 1 \
  --out "$ONE_D"
then
  SEARCH_STATUS=0
else
  SEARCH_STATUS=$?
fi

if [[ "$SEARCH_STATUS" -ne 0 && "$SEARCH_STATUS" -ne 2 ]]; then
  echo "ERROR: 1-D stage failed with exit status $SEARCH_STATUS" >&2
  exit "$SEARCH_STATUS"
fi

if [[ "$SEARCH_STATUS" -eq 2 ]]; then
  echo "V913_WEAKT_CERAMIC_1D_PARTIAL target_um=$TARGET_EXT_UM continuing_to_analysis=true" >&2
fi

"$PYTHON_BIN" scripts/analyze_v913_weakT_ceramic_1d.py \
  --case-root "$ONE_D/cases" \
  --candidate-registry "$ZERO_D/combined_promoted_registry.csv" \
  --policy-json "$SEARCH_POLICY" \
  --target-extension-um "$TARGET_EXT_UM" \
  --promote-per-class "$FINAL_PROMOTE_PER_CLASS" \
  --temperatures-K "${TEMPERATURE_ARRAY[@]}" \
  --out "$FINAL_ANALYSIS"

echo "V913_WEAKT_CERAMIC_SEARCH_COMPLETE out=$OUTROOT final_registry=$FINAL_ANALYSIS/promoted_registry.csv"
