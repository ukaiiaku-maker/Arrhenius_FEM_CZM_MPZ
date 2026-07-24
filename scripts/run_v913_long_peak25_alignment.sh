#!/usr/bin/env bash
set -euo pipefail

: "${SOURCE_RANKING:?Set SOURCE_RANKING to the prior ranked_candidates.csv}"
: "${LOADING_MAP:?Set LOADING_MAP to a mechanically calibrated long loading map}"

PYTHON_BIN="${PYTHON_BIN:-${CONDA_PREFIX:?CONDA_PREFIX is not set}/bin/python}"
BASE_PHYSICS_JSON="${BASE_PHYSICS_JSON:-mpz_v9_13_v10222_transfer_common_physics.json}"
POLICY_JSON="${POLICY_JSON:-mpz_v9_12_targeted_local_search_policy.json}"
OUTROOT="${OUTROOT:-runs/v9_13_peak25_long100um_align900K_v1}"
TARGET_EXT_UM="${TARGET_EXT_UM:-100}"
CHECKPOINTS_UM="${CHECKPOINTS_UM:-25 50 75 100}"
TEMPERATURES_K="${TEMPERATURES_K:-700 800 900 1000 1100 1200 1300}"
TARGET_PEAK_K="${TARGET_PEAK_K:-900}"
JOBS="${JOBS:-4}"
PROGRESS_INTERVAL_S="${PROGRESS_INTERVAL_S:-60}"

PREP="$OUTROOT/preparation"
SCREEN="$OUTROOT/long_screen"
ANALYSIS="$OUTROOT/alignment"
REGISTRY="$PREP/selected_peak25_registry.csv"

mkdir -p "$PREP" "$SCREEN" "$ANALYSIS"

"$PYTHON_BIN" scripts/prepare_v913_long_peak25_campaign.py \
  --source-ranking "$SOURCE_RANKING" \
  --loading-map "$LOADING_MAP" \
  --out-registry "$REGISTRY" \
  --selected-count 25 \
  --target-extension-um "$TARGET_EXT_UM"

# Invoke as a module so the repository root is always importable as `scripts`.
"$PYTHON_BIN" -u -m scripts.run_v913_autonomous_dbtt_search \
  --candidate-registry "$REGISTRY" \
  --base-physics-json "$BASE_PHYSICS_JSON" \
  --loading-map "$LOADING_MAP" \
  --policy-json "$POLICY_JSON" \
  --families \
  --per-parent 0 \
  --temperatures $TEMPERATURES_K \
  --checkpoint-um "$TARGET_EXT_UM" \
  --target-extension-um "$TARGET_EXT_UM" \
  --translation-action-exponent 0.95 \
  --max-hazard-increment 0.05 \
  --jobs "$JOBS" \
  --progress-interval-s "$PROGRESS_INTERVAL_S" \
  --promote-count 25 \
  --low-max-K 700 \
  --high-min-K 1100 \
  --peak-min-K 800 \
  --peak-max-K 1200 \
  --out "$SCREEN"

"$PYTHON_BIN" scripts/analyze_v913_long_peak_alignment.py \
  --case-root "$SCREEN/cases" \
  --candidate-registry "$REGISTRY" \
  --checkpoints-um $CHECKPOINTS_UM \
  --target-peak-temperature-K "$TARGET_PEAK_K" \
  --peak-estimator discrete \
  --stable-drift-limit-K 50 \
  --maximum-alignable-drift-K 100 \
  --minimum-post-peak-drop-MPa-sqrt-m 1 \
  --out "$ANALYSIS"

echo "V913_LONG_PEAK25_CAMPAIGN_COMPLETE out=$OUTROOT"
