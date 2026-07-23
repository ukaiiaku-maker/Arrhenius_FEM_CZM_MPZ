#!/usr/bin/env bash
set -euo pipefail

# First calibrated-1-D acquisition wave over the existing v9.12 local pool.
# The first 128 rows from each of the two peak parents form a nested,
# low-discrepancy 256-candidate Sobol prefix.  The established v9.12
# directional/peak objective and ExtraTrees acquisition are reused unchanged.

PYTHON_BIN="${PYTHON_BIN:-python}"
export PYTHONPATH="${PYTHONPATH:-.}"
REGISTRY="${REGISTRY:-candidates/v9_12_targeted_local_4096_registry.csv}"
OUT="${OUT:-runs/v9_13_autonomous_dbtt_4096_peak_wave1_v1}"
MAX_JOBS="${MAX_JOBS:-4}"
PER_PARENT="${PER_PARENT:-128}"
PROMOTE_COUNT="${PROMOTE_COUNT:-48}"
NEXT_BATCH_SIZE="${NEXT_BATCH_SIZE:-256}"
TREES="${TREES:-1200}"

PHYSICS="${PHYSICS:-mpz_v9_13_v10222_transfer_common_physics.json}"
LOADING_MAP="${LOADING_MAP:-runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_loading_map.json}"
POLICY="${POLICY:-mpz_v9_12_targeted_local_search_policy.json}"

for required in "$REGISTRY" "$PHYSICS" "$LOADING_MAP" "$POLICY"; do
  test -f "$required" || {
    echo "ERROR: missing required input: $required" >&2
    exit 1
  }
done

"$PYTHON_BIN" -m pytest -q \
  tests/test_emergent_gnd_rcurve_v913.py \
  tests/test_v913_autonomous_dbtt_search.py

"$PYTHON_BIN" -u scripts/run_v913_autonomous_dbtt_search.py \
  --candidate-registry "$REGISTRY" \
  --base-physics-json "$PHYSICS" \
  --loading-map "$LOADING_MAP" \
  --policy-json "$POLICY" \
  --families peak \
  --per-parent "$PER_PARENT" \
  --parent-offset 0 \
  --temperatures 700 800 900 1000 1100 1200 \
  --checkpoint-um 25 \
  --target-extension-um 25 \
  --translation-action-exponent 0.95 \
  --max-hazard-increment 0.25 \
  --jobs "$MAX_JOBS" \
  --promote-count "$PROMOTE_COUNT" \
  --out "$OUT"

"$PYTHON_BIN" -u scripts/train_mpz_v9_12_directional_peak_surrogate.py \
  --table "$OUT/autonomous_dbtt_training_table.csv" \
  --out-model "$OUT/autonomous_dbtt_surrogate.joblib" \
  --out-dir "$OUT/surrogate_fit" \
  --trees "$TREES" \
  --folds 5 \
  --seed 9131

"$PYTHON_BIN" -u scripts/propose_mpz_v9_12_directional_peak_batch.py \
  --model "$OUT/autonomous_dbtt_surrogate.joblib" \
  --pool-table "$OUT/candidate_pool_features.csv" \
  --pool-registry "$REGISTRY" \
  --batch-size "$NEXT_BATCH_SIZE" \
  --directional-fraction 0 \
  --peak-fraction 0.85 \
  --beta 1.5 \
  --out "$OUT/next_active_registry.csv"

echo "V913_DBTT_WAVE1_COMPLETE"
echo "ranking=$OUT/ranked_candidates.csv"
echo "promoted_registry=$OUT/promoted_registry.csv"
echo "next_active_registry=$OUT/next_active_registry.csv"
