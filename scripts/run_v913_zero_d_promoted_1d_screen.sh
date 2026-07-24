#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}}"
PYTHON_BIN="${PYTHON_BIN:-python}"
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"

SOURCE_REGISTRY="${SOURCE_REGISTRY:-runs/v9_13_zero_d_large_persistent_search_v1/promoted_registry_corrected.csv}"
BASE_PHYSICS_JSON="${BASE_PHYSICS_JSON:-mpz_v9_13_v10222_transfer_common_physics.json}"
LOADING_MAP="${LOADING_MAP:-runs/v9_13_long_map_exponential_110um_v2/v10_2_22_long_rcurve_loading_map_exponential_110um.json}"
POLICY_JSON="${POLICY_JSON:-mpz_v9_12_targeted_local_search_policy.json}"
OUTROOT="${OUTROOT:-runs/v9_13_zeroD_promoted_1d_384_50um_v1}"
SELECTED_COUNT="${SELECTED_COUNT:-384}"
TEMPERATURES_K="${TEMPERATURES_K:-700 800 900 950 1000 1050 1100 1200 1300 1400}"
TARGET_EXT_UM="${TARGET_EXT_UM:-50}"
CHECKPOINTS_UM="${CHECKPOINTS_UM:-25 50}"
JOBS="${JOBS:-4}"
PROGRESS_INTERVAL_S="${PROGRESS_INTERVAL_S:-60}"

PREP="$OUTROOT/preparation"
SCREEN="$OUTROOT/one_d_screen"
ANALYSIS="$OUTROOT/peak_analysis"
REGISTRY="$PREP/selected_${SELECTED_COUNT}_registry.csv"

for required in "$SOURCE_REGISTRY" "$BASE_PHYSICS_JSON" "$LOADING_MAP" "$POLICY_JSON"; do
  test -f "$required" || {
    echo "ERROR: missing required input: $required" >&2
    exit 1
  }
done

mkdir -p "$PREP" "$SCREEN" "$ANALYSIS"

"$PYTHON_BIN" scripts/prepare_v913_zero_d_to_1d_screen.py \
  --source-registry "$SOURCE_REGISTRY" \
  --out-registry "$REGISTRY" \
  --selected-count "$SELECTED_COUNT"

read -r -a TEMPERATURE_ARRAY <<< "$TEMPERATURES_K"
read -r -a CHECKPOINT_ARRAY <<< "$CHECKPOINTS_UM"

"$PYTHON_BIN" -u -m scripts.run_v913_autonomous_dbtt_search \
  --candidate-registry "$REGISTRY" \
  --base-physics-json "$BASE_PHYSICS_JSON" \
  --loading-map "$LOADING_MAP" \
  --policy-json "$POLICY_JSON" \
  --families \
  --per-parent 0 \
  --temperatures "${TEMPERATURE_ARRAY[@]}" \
  --checkpoint-um "$TARGET_EXT_UM" \
  --target-extension-um "$TARGET_EXT_UM" \
  --translation-action-exponent 0.95 \
  --max-hazard-increment 0.05 \
  --jobs "$JOBS" \
  --progress-interval-s "$PROGRESS_INTERVAL_S" \
  --promote-count 128 \
  --low-max-K 700 \
  --high-min-K 1200 \
  --peak-min-K 850 \
  --peak-max-K 1100 \
  --direction-threshold 5 \
  --peak-threshold 5 \
  --out "$SCREEN"

"$PYTHON_BIN" scripts/analyze_v913_long_peak_alignment.py \
  --case-root "$SCREEN/cases" \
  --candidate-registry "$REGISTRY" \
  --checkpoints-um "${CHECKPOINT_ARRAY[@]}" \
  --target-peak-temperature-K 900 \
  --peak-estimator discrete \
  --stable-drift-limit-K 50 \
  --maximum-alignable-drift-K 100 \
  --minimum-post-peak-drop-MPa-sqrt-m 5 \
  --refinement-step-K 25 \
  --refinement-half-width-K 100 \
  --out "$ANALYSIS"

echo "V913_ZERO_D_PROMOTED_1D_SCREEN_COMPLETE out=$OUTROOT registry=$REGISTRY"
