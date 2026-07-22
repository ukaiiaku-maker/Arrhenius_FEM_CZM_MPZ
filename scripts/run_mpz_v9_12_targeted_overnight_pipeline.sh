#!/usr/bin/env bash
set -euo pipefail

# v9.12 targeted 0-D -> 1-D overnight acquisition pipeline.
# Expected to be run from the repository root in arrhenius-fem-czm-v912-gnd.

REGISTRY_20K="${REGISTRY_20K:-candidates/v9_12_state_focused_overnight_r0_20000.csv}"
ROOT_0D_20K="${ROOT_0D_20K:-runs/v9_12_state_focused_0d_overnight_r0_20000_v1}"
REGISTRY_384="${REGISTRY_384:-candidates/v9_12_1d_balanced_promotion_384_registry.csv}"
ROOT_1D_384="${ROOT_1D_384:-runs/v9_12_1d_balanced_promotion_384_coupled_c2_dt0p1_v1}"

WORK="${WORK:-ml/v9_12_targeted_directional_peak_overnight_v1}"
LOCAL_REGISTRY="${LOCAL_REGISTRY:-candidates/v9_12_targeted_local_4096_registry.csv}"
LOCAL_0D_ROOT="${LOCAL_0D_ROOT:-runs/v9_12_targeted_local_4096_0d_v1}"
LOCAL_BATCH="${LOCAL_BATCH:-candidates/v9_12_targeted_local_active_1152_registry.csv}"
GLOBAL_BATCH="${GLOBAL_BATCH:-candidates/v9_12_targeted_global_active_384_registry.csv}"
COMBINED_REGISTRY="${COMBINED_REGISTRY:-candidates/v9_12_targeted_overnight_1536_registry.csv}"
SMOKE_REGISTRY="${SMOKE_REGISTRY:-candidates/v9_12_targeted_overnight_smoke_8_registry.csv}"
SMOKE_ROOT="${SMOKE_ROOT:-runs/v9_12_targeted_overnight_smoke_8_v1}"
OUT_1D="${OUT_1D:-runs/v9_12_targeted_directional_peak_1d_1536_coupled_c2_dt0p1_v1}"

LOCAL_PER_SEED="${LOCAL_PER_SEED:-512}"
LOCAL_BATCH_SIZE="${LOCAL_BATCH_SIZE:-1152}"
GLOBAL_BATCH_SIZE="${GLOBAL_BATCH_SIZE:-384}"
EXPECTED_COMBINED="${EXPECTED_COMBINED:-1536}"

TEMPERATURES=(300 400 500 600 700 800 900 1000 1100 1200)
PROTOCOL="mpz_v9_12_protocol_example.csv"
PHYSICS="mpz_v9_12_emergent_gnd_common_physics.json"
BOUNDS="mpz_v9_12_emergent_gnd_search_bounds_state_focused.json"
POLICY="mpz_v9_12_targeted_local_search_policy.json"

mkdir -p "$WORK" candidates runs
PIPELINE_LOG="$WORK/pipeline.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

echo "TARGETED_OVERNIGHT_START $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "git_head=$(git rev-parse HEAD)"

for required in \
  "$REGISTRY_20K" "$ROOT_0D_20K" \
  "$REGISTRY_384" "$ROOT_1D_384" \
  "$PROTOCOL" "$PHYSICS" "$BOUNDS" "$POLICY"; do
  test -e "$required" || {
    echo "ERROR: required input missing: $required" >&2
    exit 1
  }
done

python -m py_compile \
  scripts/generate_mpz_v9_12_targeted_local_pool.py \
  scripts/build_mpz_v9_12_0d_to_1d_transfer_table.py \
  scripts/augment_mpz_v9_12_directional_peak_targets.py \
  scripts/train_mpz_v9_12_directional_peak_surrogate.py \
  scripts/propose_mpz_v9_12_directional_peak_batch.py \
  scripts/run_mpz_v9_12_emergent_gnd_screen_resilient.py

EMPTY_1D="$WORK/empty_1d"
mkdir -p "$EMPTY_1D"

TRAIN_BASE="$WORK/training_384_base.csv"
TRAIN_AUG="$WORK/training_384_directional_peak.csv"
MODEL="$WORK/directional_peak_surrogate.joblib"
FIT_DIR="$WORK/fit"

python -u scripts/build_mpz_v9_12_0d_to_1d_transfer_table.py \
  --candidate-registry "$REGISTRY_384" \
  --zero-d-root "$ROOT_0D_20K" \
  --one-d-root "$ROOT_1D_384" \
  --bounds-json "$BOUNDS" \
  --out "$TRAIN_BASE"

python -u scripts/augment_mpz_v9_12_directional_peak_targets.py \
  --input "$TRAIN_BASE" \
  --out "$TRAIN_AUG" \
  --direction-threshold 5 \
  --peak-threshold 1

python -u scripts/train_mpz_v9_12_directional_peak_surrogate.py \
  --table "$TRAIN_AUG" \
  --out-model "$MODEL" \
  --out-dir "$FIT_DIR" \
  --trees 1200 \
  --folds 5 \
  --seed 9122

python -u scripts/generate_mpz_v9_12_targeted_local_pool.py \
  --base-registry "$REGISTRY_20K" \
  --policy-json "$POLICY" \
  --per-seed "$LOCAL_PER_SEED" \
  --seed 19120 \
  --prefix v912_targeted_local \
  --out "$LOCAL_REGISTRY"

run_resilient () {
  local stage="$1"
  local registry="$2"
  local out="$3"
  local progress_every="$4"
  local compact="$5"
  local -a resume_args=()
  local -a compact_args=()

  if test -f "$out/resilient_records_checkpoint.csv"; then
    resume_args+=(--resume)
  fi
  if test "$compact" = "1"; then
    compact_args+=(--compact-output)
  fi

  caffeinate -dimsu env \
  PYTHONUNBUFFERED=1 \
  MPZ_V912_COUPLED_OPERATOR_SUBSTEPS=2 \
  MPZ_V912_MAX_FEEDBACK_SUBSTEP_S=0.1 \
  python -u scripts/run_mpz_v9_12_emergent_gnd_screen_resilient.py \
    "${resume_args[@]}" \
    --stage "$stage" \
    --candidate-registry "$registry" \
    --protocol-csv "$PROTOCOL" \
    --physics-json "$PHYSICS" \
    --temperatures "${TEMPERATURES[@]}" \
    --window-um 10 30 \
    --min-amplitude 50 \
    --target-localization 0.50 \
    --max-width-K 200 \
    "${compact_args[@]}" \
    --quiet-inner \
    --progress-every "$progress_every" \
    --out "$out"
}

run_resilient 0d "$LOCAL_REGISTRY" "$LOCAL_0D_ROOT" 64 1

LOCAL_POOL_BASE="$WORK/local_pool_base.csv"
LOCAL_POOL_AUG="$WORK/local_pool_directional_peak.csv"
GLOBAL_POOL_BASE="$WORK/global_pool_base.csv"
GLOBAL_POOL_AUG="$WORK/global_pool_directional_peak.csv"

python -u scripts/build_mpz_v9_12_0d_to_1d_transfer_table.py \
  --candidate-registry "$LOCAL_REGISTRY" \
  --zero-d-root "$LOCAL_0D_ROOT" \
  --one-d-root "$EMPTY_1D" \
  --bounds-json "$BOUNDS" \
  --out "$LOCAL_POOL_BASE"

python -u scripts/augment_mpz_v9_12_directional_peak_targets.py \
  --input "$LOCAL_POOL_BASE" \
  --out "$LOCAL_POOL_AUG"

python -u scripts/build_mpz_v9_12_0d_to_1d_transfer_table.py \
  --candidate-registry "$REGISTRY_20K" \
  --zero-d-root "$ROOT_0D_20K" \
  --one-d-root "$EMPTY_1D" \
  --bounds-json "$BOUNDS" \
  --out "$GLOBAL_POOL_BASE"

python -u scripts/augment_mpz_v9_12_directional_peak_targets.py \
  --input "$GLOBAL_POOL_BASE" \
  --out "$GLOBAL_POOL_AUG"

python -u scripts/propose_mpz_v9_12_directional_peak_batch.py \
  --model "$MODEL" \
  --pool-table "$LOCAL_POOL_AUG" \
  --pool-registry "$LOCAL_REGISTRY" \
  --batch-size "$LOCAL_BATCH_SIZE" \
  --directional-fraction 0.50 \
  --peak-fraction 0.375 \
  --beta 1.5 \
  --out "$LOCAL_BATCH"

python -u scripts/propose_mpz_v9_12_directional_peak_batch.py \
  --model "$MODEL" \
  --pool-table "$GLOBAL_POOL_AUG" \
  --pool-registry "$REGISTRY_20K" \
  --exclude-registry "$REGISTRY_384" \
  --batch-size "$GLOBAL_BATCH_SIZE" \
  --directional-fraction 0.5833333333333334 \
  --peak-fraction 0.25 \
  --beta 2.0 \
  --out "$GLOBAL_BATCH"

python - "$LOCAL_BATCH" "$GLOBAL_BATCH" "$COMBINED_REGISTRY" "$EXPECTED_COMBINED" <<'PY'
from pathlib import Path
import sys
import pandas as pd

local_path, global_path, out_path, expected_text = sys.argv[1:]
local = pd.read_csv(local_path, low_memory=False)
global_ = pd.read_csv(global_path, low_memory=False)
local["overnight_source_pool"] = "targeted_local"
global_["overnight_source_pool"] = "global_20k"
combined = pd.concat([local, global_], ignore_index=True, sort=False)
expected = int(expected_text)
if len(combined) != expected:
    raise RuntimeError(f"expected {expected} rows, found {len(combined)}")
if combined["candidate_id"].nunique() != expected:
    raise RuntimeError("combined registry contains duplicate candidate IDs")
Path(out_path).parent.mkdir(parents=True, exist_ok=True)
combined.to_csv(out_path, index=False)
print("COMBINED_REGISTRY", len(combined), out_path)
print(combined.groupby(["overnight_source_pool", "acquisition_role"]).size())
PY

python - "$COMBINED_REGISTRY" "$SMOKE_REGISTRY" <<'PY'
from pathlib import Path
import sys
import pandas as pd

source_path, out_path = sys.argv[1:]
df = pd.read_csv(source_path, low_memory=False)
parts = []
for role, count in (("directional", 3), ("peak", 3), ("exploration", 2)):
    subset = df[df["acquisition_role"].eq(role)].head(count)
    if len(subset) != count:
        raise RuntimeError(f"insufficient {role} candidates for smoke")
    parts.append(subset)
smoke = pd.concat(parts, ignore_index=True, sort=False)
Path(out_path).parent.mkdir(parents=True, exist_ok=True)
smoke.to_csv(out_path, index=False)
print("SMOKE_REGISTRY", len(smoke), out_path)
PY

run_resilient 1d "$SMOKE_REGISTRY" "$SMOKE_ROOT" 1 0

python - "$SMOKE_ROOT/resilient_progress.json" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
print(json.dumps(payload, indent=2, sort_keys=True))
if payload.get("complete_count") != payload.get("total_candidates"):
    raise RuntimeError("targeted 1-D smoke did not complete cleanly")
PY

run_resilient 1d "$COMBINED_REGISTRY" "$OUT_1D" 16 1

echo "TARGETED_OVERNIGHT_COMPLETE $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "output=$OUT_1D"
