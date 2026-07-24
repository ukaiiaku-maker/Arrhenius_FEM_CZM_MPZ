#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}}"
PYTHON_BIN="${PYTHON_BIN:-python}"
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"

ANCHOR_REGISTRY="${ANCHOR_REGISTRY:-candidates/v9_12_targeted_local_4096_registry.csv}"
PHYSICS="${PHYSICS:-mpz_v9_13_v10222_transfer_common_physics.json}"
LOADING_MAP="${LOADING_MAP:-runs/v9_13_long_map_exponential_110um_v2/v10_2_22_long_rcurve_loading_map_exponential_110um.json}"
POLICY="${POLICY:-mpz_v9_13_zero_d_large_search_policy.json}"
OUT="${OUT:-runs/v9_13_zero_d_large_persistent_search_v1}"

SAMPLES="${SAMPLES:-262144}"
PROXY_BATCH_SIZE="${PROXY_BATCH_SIZE:-4096}"
EXACT_COUNT="${EXACT_COUNT:-4096}"
PROMOTE_COUNT="${PROMOTE_COUNT:-512}"
MAX_JOBS="${MAX_JOBS:-4}"
SEED="${SEED:-913100}"
TEMPERATURES_K="${TEMPERATURES_K:-700 800 900 950 1000 1050 1100 1200 1300 1400}"
PROXY_EXTENSION_UM="${PROXY_EXTENSION_UM:-50}"
EXACT_EXTENSION_UM="${EXACT_EXTENSION_UM:-50}"
CHECKPOINT_UM="${CHECKPOINT_UM:-50}"
LOAD_INCREMENT_FACTOR="${LOAD_INCREMENT_FACTOR:-2}"
PROGRESS_INTERVAL_S="${PROGRESS_INTERVAL_S:-60}"

if [[ ! -f "$ANCHOR_REGISTRY" ]]; then
  "$PYTHON_BIN" scripts/materialize_v913_candidate_registry.py \
    --out "$ANCHOR_REGISTRY"
fi

for required in "$ANCHOR_REGISTRY" "$PHYSICS" "$LOADING_MAP" "$POLICY"; do
  test -f "$required" || {
    echo "ERROR: missing required input: $required" >&2
    exit 1
  }
done

"$PYTHON_BIN" -m pytest -q \
  tests/test_zero_d_persistent_v913.py \
  tests/test_v913_zero_d_large_search.py

mkdir -p "$OUT"

read -r -a TEMPERATURE_ARRAY <<< "$TEMPERATURES_K"

"$PYTHON_BIN" -u scripts/run_v913_zero_d_large_search.py \
  --anchor-registry "$ANCHOR_REGISTRY" \
  --base-physics-json "$PHYSICS" \
  --loading-map "$LOADING_MAP" \
  --policy-json "$POLICY" \
  --out "$OUT" \
  --samples "$SAMPLES" \
  --proxy-batch-size "$PROXY_BATCH_SIZE" \
  --exact-count "$EXACT_COUNT" \
  --promote-count "$PROMOTE_COUNT" \
  --jobs "$MAX_JOBS" \
  --seed "$SEED" \
  --temperatures-K "${TEMPERATURE_ARRAY[@]}" \
  --proxy-extension-um "$PROXY_EXTENSION_UM" \
  --exact-extension-um "$EXACT_EXTENSION_UM" \
  --checkpoint-um "$CHECKPOINT_UM" \
  --load-increment-factor "$LOAD_INCREMENT_FACTOR" \
  --progress-interval-s "$PROGRESS_INTERVAL_S"
