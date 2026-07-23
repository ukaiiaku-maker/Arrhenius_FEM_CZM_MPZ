#!/usr/bin/env bash
set -euo pipefail

# First calibrated-1-D acquisition wave over the existing v9.12 local pool.
# The first 128 rows from each of the two peak parents form a nested,
# low-discrepancy 256-candidate Sobol prefix.  The established v9.12
# directional/peak objective and ExtraTrees acquisition are reused unchanged.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python}"
# The checked-out branch must take precedence over any older editable install.
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"
DEFAULT_REGISTRY="candidates/v9_12_targeted_local_4096_registry.csv"
REGISTRY="${REGISTRY:-$DEFAULT_REGISTRY}"
OUT="${OUT:-runs/v9_13_autonomous_dbtt_4096_peak_wave1_v1}"
MAX_JOBS="${MAX_JOBS:-4}"
PER_PARENT="${PER_PARENT:-128}"
PROMOTE_COUNT="${PROMOTE_COUNT:-48}"
NEXT_BATCH_SIZE="${NEXT_BATCH_SIZE:-256}"
TREES="${TREES:-1200}"
PROGRESS_INTERVAL_S="${PROGRESS_INTERVAL_S:-60}"

PHYSICS="${PHYSICS:-mpz_v9_13_v10222_transfer_common_physics.json}"
LOADING_MAP="${LOADING_MAP:-runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_loading_map.json}"
POLICY="${POLICY:-mpz_v9_12_targeted_local_search_policy.json}"

mkdir -p "$OUT"
WAVE_STATUS="$OUT/wave1_status.json"
WAVE_STARTED_AT="$("$PYTHON_BIN" - <<'PY'
from datetime import datetime, timezone
print(datetime.now(timezone.utc).isoformat())
PY
)"
CURRENT_STAGE="initialization"
WAVE_COMPLETE=0

write_wave_status() {
  local state="$1"
  local stage_status="$2"
  local exit_code="${3:-}"
  "$PYTHON_BIN" - \
    "$WAVE_STATUS" \
    "$CURRENT_STAGE" \
    "$state" \
    "$stage_status" \
    "$WAVE_STARTED_AT" \
    "$exit_code" <<'PY'
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

path = Path(sys.argv[1])
payload = {
    "schema": "v9.13_dbtt_wave1_status_v1",
    "stage": sys.argv[2],
    "state": sys.argv[3],
    "stage_status": sys.argv[4],
    "started_at_utc": sys.argv[5],
    "updated_at_utc": datetime.now(timezone.utc).isoformat(),
    "exit_code": int(sys.argv[6]) if sys.argv[6] else None,
}
temporary = path.with_suffix(path.suffix + ".tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
temporary.replace(path)
PY
}

stage_start() {
  CURRENT_STAGE="$1"
  write_wave_status "running" "started"
  echo "V913_DBTT_STAGE stage=$CURRENT_STAGE status=start"
}

stage_complete() {
  write_wave_status "running" "complete"
  echo "V913_DBTT_STAGE stage=$CURRENT_STAGE status=complete"
}

on_exit() {
  local exit_code=$?
  if [[ "$WAVE_COMPLETE" -ne 1 ]]; then
    set +e
    write_wave_status "failed" "failed" "$exit_code"
    echo \
      "V913_DBTT_WAVE1_FAILED stage=$CURRENT_STAGE exit_code=$exit_code" \
      >&2
  fi
}
trap on_exit EXIT

stage_start "candidate_registry"
if [[ "$REGISTRY" == "$DEFAULT_REGISTRY" ]]; then
  "$PYTHON_BIN" -u scripts/materialize_v913_candidate_registry.py \
    --out "$REGISTRY"
fi
stage_complete

stage_start "required_inputs"
for required in "$REGISTRY" "$PHYSICS" "$LOADING_MAP" "$POLICY"; do
  test -f "$required" || {
    echo "ERROR: missing required input: $required" >&2
    exit 1
  }
done
stage_complete

stage_start "import_validation"
"$PYTHON_BIN" - "$REPO_ROOT" <<'PY'
from pathlib import Path
import sys

import arrhenius_fracture

expected = Path(sys.argv[1]).resolve()
loaded = Path(arrhenius_fracture.__file__).resolve()
if expected not in loaded.parents:
    raise SystemExit(
        f"ERROR: arrhenius_fracture loaded from {loaded}, expected {expected}"
    )
print(f"V913_DBTT_IMPORT_OK module={loaded}")
PY
stage_complete

stage_start "focused_tests"
"$PYTHON_BIN" -m pytest -q \
  tests/test_emergent_gnd_rcurve_v913.py \
  tests/test_v913_autonomous_dbtt_search.py \
  tests/test_v913_candidate_registry.py
stage_complete

stage_start "integration_sentinels"
"$PYTHON_BIN" -u scripts/verify_v913_autonomous_dbtt_integration.py \
  --candidate-registry "$REGISTRY" \
  --base-physics-json "$PHYSICS" \
  --loading-map "$LOADING_MAP" \
  --run-sentinels
stage_complete

stage_start "candidate_search"
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
  --max-hazard-increment 0.05 \
  --jobs "$MAX_JOBS" \
  --progress-interval-s "$PROGRESS_INTERVAL_S" \
  --promote-count "$PROMOTE_COUNT" \
  --out "$OUT"
stage_complete

stage_start "surrogate_fit"
"$PYTHON_BIN" -u scripts/train_mpz_v9_12_directional_peak_surrogate.py \
  --table "$OUT/autonomous_dbtt_training_table.csv" \
  --out-model "$OUT/autonomous_dbtt_surrogate.joblib" \
  --out-dir "$OUT/surrogate_fit" \
  --trees "$TREES" \
  --folds 5 \
  --seed 9131
stage_complete

stage_start "active_batch_proposal"
"$PYTHON_BIN" -u scripts/propose_mpz_v9_12_directional_peak_batch.py \
  --model "$OUT/autonomous_dbtt_surrogate.joblib" \
  --pool-table "$OUT/candidate_pool_features.csv" \
  --pool-registry "$REGISTRY" \
  --batch-size "$NEXT_BATCH_SIZE" \
  --directional-fraction 0 \
  --peak-fraction 0.85 \
  --beta 1.5 \
  --out "$OUT/next_active_registry.csv"
stage_complete

CURRENT_STAGE="complete"
write_wave_status "complete" "complete" "0"
WAVE_COMPLETE=1
echo "V913_DBTT_WAVE1_COMPLETE"
echo "ranking=$OUT/ranked_candidates.csv"
echo "promoted_registry=$OUT/promoted_registry.csv"
echo "next_active_registry=$OUT/next_active_registry.csv"
