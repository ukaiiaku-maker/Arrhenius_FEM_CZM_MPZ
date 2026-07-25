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
ZERO_D_LOADING_MAP="${ZERO_D_LOADING_MAP:-runs/v9_13_long_map_exponential_110um_v2/v10_2_22_long_rcurve_loading_map_exponential_110um.json}"
LONG_LOADING_MAP="${LONG_LOADING_MAP:-runs/v9_13_long_map_exponential_1000um_v1/v10_2_22_long_rcurve_loading_map_exponential_1000um.json}"
SEARCH_POLICY="${SEARCH_POLICY:-mpz_v9_13_zero_d_weakT_ceramic_search_policy.json}"
LEGACY_POLICY="${LEGACY_POLICY:-mpz_v9_12_targeted_local_search_policy.json}"
OUTROOT="${OUTROOT:-runs/v9_13_weakT_ceramic_search_v1}"
TEMPERATURES_K="${TEMPERATURES_K:-300 400 500 600 700 800 900 950 1000 1050 1100 1150 1200 1250 1300}"
SAMPLES="${SAMPLES:-262144}"
ZERO_D_EXACT_PER_CLASS="${ZERO_D_EXACT_PER_CLASS:-512}"
ZERO_D_PROMOTE_PER_CLASS="${ZERO_D_PROMOTE_PER_CLASS:-48}"
STAGE1_PROMOTE_PER_CLASS="${STAGE1_PROMOTE_PER_CLASS:-10}"
FINAL_PROMOTE_PER_CLASS="${FINAL_PROMOTE_PER_CLASS:-5}"
MAX_JOBS="${MAX_JOBS:-4}"
PROGRESS_INTERVAL_S="${PROGRESS_INTERVAL_S:-60}"

ZERO_D="$OUTROOT/zero_d"
STAGE1="$OUTROOT/one_d_250um"
STAGE1_ANALYSIS="$OUTROOT/analysis_250um"
STAGE2="$OUTROOT/one_d_1000um"
FINAL_ANALYSIS="$OUTROOT/analysis_1000um"

for required in \
  "$ANCHOR_REGISTRY" \
  "$BASE_PHYSICS_JSON" \
  "$ZERO_D_LOADING_MAP" \
  "$LONG_LOADING_MAP" \
  "$SEARCH_POLICY" \
  "$LEGACY_POLICY"
do
  test -f "$required" || {
    echo "ERROR: missing required input: $required" >&2
    exit 1
  }
done

mkdir -p "$ZERO_D" "$STAGE1" "$STAGE1_ANALYSIS" "$STAGE2" "$FINAL_ANALYSIS"
read -r -a TEMPERATURE_ARRAY <<< "$TEMPERATURES_K"

"$PYTHON_BIN" - "$ZERO_D_LOADING_MAP" "$LONG_LOADING_MAP" <<'PY'
import json
import sys
from pathlib import Path

for raw, required in ((sys.argv[1], 100.0), (sys.argv[2], 1000.0)):
    path = Path(raw)
    payload = json.loads(path.read_text())
    coverage = sum(float(value) for value in payload["projected_advances_m"]) * 1.0e6
    if coverage + 1.0e-9 < required:
        raise SystemExit(
            f"ERROR: loading map {path} covers {coverage:.6g} um; required {required:g} um"
        )
    print(f"LOADING_MAP_OK path={path} coverage_um={coverage:.9g}")
PY

"$PYTHON_BIN" -m pytest -q \
  tests/test_v913_weakT_ceramic_objectives.py \
  tests/test_v913_weakT_ceramic_pipeline.py

"$PYTHON_BIN" -u scripts/run_v913_zero_d_weakT_ceramic_search.py \
  --anchor-registry "$ANCHOR_REGISTRY" \
  --base-physics-json "$BASE_PHYSICS_JSON" \
  --loading-map "$ZERO_D_LOADING_MAP" \
  --policy-json "$SEARCH_POLICY" \
  --samples "$SAMPLES" \
  --exact-per-class "$ZERO_D_EXACT_PER_CLASS" \
  --promote-per-class "$ZERO_D_PROMOTE_PER_CLASS" \
  --jobs "$MAX_JOBS" \
  --temperatures-K "${TEMPERATURE_ARRAY[@]}" \
  --progress-interval-s "$PROGRESS_INTERVAL_S" \
  --out "$ZERO_D"

run_one_d() {
  local registry="$1"
  local target_um="$2"
  local out="$3"
  local status=0
  if "$PYTHON_BIN" -u -m scripts.run_v913_autonomous_dbtt_search \
    --candidate-registry "$registry" \
    --base-physics-json "$BASE_PHYSICS_JSON" \
    --loading-map "$LONG_LOADING_MAP" \
    --policy-json "$LEGACY_POLICY" \
    --families \
    --per-parent 0 \
    --temperatures "${TEMPERATURE_ARRAY[@]}" \
    --checkpoint-um "$target_um" \
    --target-extension-um "$target_um" \
    --translation-action-exponent 0.95 \
    --max-hazard-increment 0.05 \
    --jobs "$MAX_JOBS" \
    --progress-interval-s "$PROGRESS_INTERVAL_S" \
    --promote-count 1 \
    --out "$out"
  then
    status=0
  else
    status=$?
  fi
  if [[ "$status" -ne 0 && "$status" -ne 2 ]]; then
    echo "ERROR: 1-D stage failed with exit status $status" >&2
    exit "$status"
  fi
  if [[ "$status" -eq 2 ]]; then
    echo "V913_WEAKT_CERAMIC_1D_PARTIAL target_um=$target_um continuing_to_analysis=true" >&2
  fi
}

run_one_d "$ZERO_D/combined_promoted_registry.csv" 250 "$STAGE1"

"$PYTHON_BIN" scripts/analyze_v913_weakT_ceramic_1d.py \
  --case-root "$STAGE1/cases" \
  --candidate-registry "$ZERO_D/combined_promoted_registry.csv" \
  --policy-json "$SEARCH_POLICY" \
  --target-extension-um 250 \
  --promote-per-class "$STAGE1_PROMOTE_PER_CLASS" \
  --temperatures-K "${TEMPERATURE_ARRAY[@]}" \
  --out "$STAGE1_ANALYSIS"

run_one_d "$STAGE1_ANALYSIS/promoted_registry.csv" 1000 "$STAGE2"

"$PYTHON_BIN" scripts/analyze_v913_weakT_ceramic_1d.py \
  --case-root "$STAGE2/cases" \
  --candidate-registry "$STAGE1_ANALYSIS/promoted_registry.csv" \
  --policy-json "$SEARCH_POLICY" \
  --target-extension-um 1000 \
  --promote-per-class "$FINAL_PROMOTE_PER_CLASS" \
  --temperatures-K "${TEMPERATURE_ARRAY[@]}" \
  --out "$FINAL_ANALYSIS"

echo "V913_WEAKT_CERAMIC_SEARCH_COMPLETE out=$OUTROOT final_registry=$FINAL_ANALYSIS/promoted_registry.csv"
