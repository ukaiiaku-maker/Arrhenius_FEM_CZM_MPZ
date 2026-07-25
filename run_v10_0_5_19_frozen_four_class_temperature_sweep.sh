#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_19_frozen4_anisotropic_300_1200K_1000um_v1}
FAMILY_JSON=${FAMILY_JSON:-$ROOT/runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json}
FAMILY_JSON_SHA256_REQUIRED=${FAMILY_JSON_SHA256_REQUIRED:-a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde}

TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 850 900 950 1000 1050 1100 1150 1200"}
TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
STEPS=${STEPS:-1000000}
DU=${DU:-2e-5}
DT=${DT:-840}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
RESTART_INCOMPLETE=${RESTART_INCOMPLETE:-1}
PRINT_EVERY=${PRINT_EVERY:-100}
BASE_HAZARD_SEED=${BASE_HAZARD_SEED:-1720}
SEED_NAMESPACE=${SEED_NAMESPACE:-v10.0.5.19-frozen-four-class-production-v1}
EVENT_MIN_FACTOR=${EVENT_MIN_FACTOR:-0.5}
EVENT_MAX_FACTOR=${EVENT_MAX_FACTOR:-4.0}
QUALITY_SUBDIVISION_COUNTS=${QUALITY_SUBDIVISION_COUNTS:-2,3,4,6,8}
QUALITY_PARTITION_MAX_DEPTH=${QUALITY_PARTITION_MAX_DEPTH:-9}
QUALITY_PARTITION_MIN_SEGMENT_M=${QUALITY_PARTITION_MIN_SEGMENT_M:-2e-8}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-20}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-1}
CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-30}
MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}

OPTIONS=(
  v913_paper_peak01_0242980_persistent_sites
  v913_paper_dbtt01_0202500_persistent_sites
  v913_paper_weakT01_0129902_persistent_sites
  v913_paper_ceramic01_0077080_persistent_sites
)

for VALUE_NAME in SKIP_FINISHED RESTART_INCOMPLETE GENERATE_SOLVER_PLOTS; do
  VALUE=${!VALUE_NAME}
  case "$VALUE" in
    0|1) ;;
    *) echo "ERROR: $VALUE_NAME must be 0 or 1" >&2; exit 1 ;;
  esac
done

if [[ ! -f "$FAMILY_JSON" ]]; then
  echo "ERROR: missing local signed-kernel family: $FAMILY_JSON" >&2
  exit 1
fi

FAMILY_JSON=$(cd "$(dirname "$FAMILY_JSON")" && pwd)/$(basename "$FAMILY_JSON")
FAMILY_JSON_SHA256=$(python - "$FAMILY_JSON" <<'PY'
from pathlib import Path
import hashlib
import sys
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)

if [[ "$FAMILY_JSON_SHA256" != "$FAMILY_JSON_SHA256_REQUIRED" ]]; then
  echo "ERROR: signed-kernel family SHA-256 mismatch" >&2
  echo "Expected: $FAMILY_JSON_SHA256_REQUIRED" >&2
  echo "Observed: $FAMILY_JSON_SHA256" >&2
  exit 1
fi

mkdir -p "$CAMPAIGN_ROOT"
echo "$$" > "$CAMPAIGN_ROOT/runner.pid"

python "$ROOT/scripts/validate_v100517_frozen_four_class.py" \
  --out "$CAMPAIGN_ROOT/frozen_parameter_catalog.json"

CASE_MATRIX="$CAMPAIGN_ROOT/case_matrix.tsv"
SEED_MAP="$CAMPAIGN_ROOT/parameter_temperature_seed_map.csv"
OPTIONS_TEXT=$(printf '%s\n' "${OPTIONS[@]}")
OPTIONS_TEXT="$OPTIONS_TEXT" TEMPERATURES="$TEMPERATURES" \
BASE_HAZARD_SEED="$BASE_HAZARD_SEED" SEED_NAMESPACE="$SEED_NAMESPACE" \
CASE_MATRIX="$CASE_MATRIX" SEED_MAP="$SEED_MAP" python - <<'PY'
import csv
import hashlib
import os
from pathlib import Path

options = [line.strip() for line in os.environ["OPTIONS_TEXT"].splitlines() if line.strip()]
temperatures = [int(value) for value in os.environ["TEMPERATURES"].split()]
base = int(os.environ["BASE_HAZARD_SEED"])
namespace = os.environ["SEED_NAMESPACE"]
rows = []
seen = set()
for option in options:
    for temperature in temperatures:
        token = f"{namespace}|{base}|{option}|{temperature}".encode()
        seed = int.from_bytes(hashlib.sha256(token).digest()[:4], "big")
        if seed in seen:
            raise SystemExit(f"seed collision for {option} at {temperature} K")
        seen.add(seed)
        rows.append(("v10.0.5.17-frozen-four-class", option, temperature, seed))
Path(os.environ["CASE_MATRIX"]).write_text(
    "".join(f"{a}\t{b}\t{c}\t{d}\n" for a, b, c, d in rows)
)
with Path(os.environ["SEED_MAP"]).open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["parameter_entry", "parameter_option", "temperature_K", "hazard_seed"])
    writer.writerows(rows)
PY

EXPECTED_CASES=$(wc -l < "$CASE_MATRIX" | tr -d ' ')
if [[ "$EXPECTED_CASES" -ne 56 ]]; then
  echo "ERROR: expected 56 cases, found $EXPECTED_CASES" >&2
  exit 1
fi

cat > "$CAMPAIGN_ROOT/campaign_configuration.txt" <<EOF
release=10.0.5.19-frozen-four-class-adaptive-quality-partition
entry=arrhenius_fracture.mode_i_first_passage_v10_0_5_19_four_class_standalone
parameter_source=frozen_FEM_registry
phase_field_solver_dependency=false
phase_field_result_comparison_enabled=false
parameter_options=${OPTIONS[*]}
temperatures_K=$TEMPERATURES
expected_cases=$EXPECTED_CASES
target_extension_um=$TARGET_EXT_UM
steps=$STEPS
dU_m=$DU
dt_s=$DT
base_hazard_seed=$BASE_HAZARD_SEED
seed_namespace=$SEED_NAMESPACE
hazard_seed_formula=sha256(seed_namespace|base|parameter_option|temperature_K)_first_32_bits
cleavage_hazard_mode=minus_log_uniform_exponential_integrated_action
cleavage_event_length_mode=threshold_scaled_same_Xi_mean_preserving
event_min_factor=$EVENT_MIN_FACTOR
event_max_factor=$EVENT_MAX_FACTOR
adaptive_czm_uniform_subdivision_fast_path=$QUALITY_SUBDIVISION_COUNTS
adaptive_czm_recursive_partition_max_depth=$QUALITY_PARTITION_MAX_DEPTH
adaptive_czm_recursive_partition_min_segment_m=$QUALITY_PARTITION_MIN_SEGMENT_M
triangle_quality_floor_relaxed=false
child_area_ratio_floor_relaxed=false
snapshot_spacing_um=$SNAPSHOT_BY_EXT_UM
save_snapshots=$SAVE_SNAPSHOTS
snapshot_cols=$SNAPSHOT_COLS
generate_solver_plots=$GENERATE_SOLVER_PLOTS
crystal_anisotropic_elasticity=true
crystal_direction_competition=true
crystal_theta_deg=$CRYSTAL_THETA_DEG
crystal_plane_gate=global_ligament_forward
minimum_global_forward=$MIN_GLOBAL_FORWARD
branching=false
maximum_fronts=1
kernel_family=$FAMILY_JSON
kernel_family_sha256=$FAMILY_JSON_SHA256
max_jobs=$MAX_JOBS
skip_finished=$SKIP_FINISHED
restart_incomplete=$RESTART_INCOMPLETE
EOF

run_case() {
  set -euo pipefail
  local ENTRY="$1"
  local OPTION="$2"
  local T="$3"
  local SEED="$4"
  local TAG
  printf -v TAG "%04d" "$T"
  local OUT="$CAMPAIGN_ROOT/$OPTION/T${TAG}K"
  local MANIFEST="$OUT/persistent_site_production_manifest_v10_0_5_17.json"
  local LOG="$OUT/console.log"

  if [[ "$SKIP_FINISHED" == "1" && -f "$MANIFEST" ]] &&
     grep -q '"run_completed_without_exception": true' "$MANIFEST"; then
    echo "[SKIP] option=$OPTION T=${T}K"
    return 0
  fi

  if [[ -d "$OUT" ]]; then
    if [[ "$RESTART_INCOMPLETE" == "1" ]]; then
      echo "[RESTART] option=$OPTION T=${T}K"
      rm -rf "$OUT"
    else
      echo "ERROR: incomplete output exists and RESTART_INCOMPLETE=0: $OUT" >&2
      return 1
    fi
  fi

  mkdir -p "$OUT"
  echo "[START] option=$OPTION entry=$ENTRY T=${T}K seed=$SEED"

  local -a CMD
  CMD=(
    python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_19_four_class_standalone
    --parameter-entry "$ENTRY"
    --parameter-option "$OPTION"
    --signed-kernel-family "$FAMILY_JSON"
    --tip-refinement-radius-um 330
    --selected-cluster-J-outer-um 240
    --local-J-outer-um 100
    --mode 2d
    --bulk-plasticity-mode tip_only
    --temperatures "$T"
    --steps "$STEPS"
    --nx 36 --ny 72
    --tip-h-fine 2.5e-6 --tip-ratio 1.15
    --dU "$DU" --dt "$DT"
    --n-stagger 1
    --print-every "$PRINT_EVERY"
    --adaptive-events
    --adaptive-event-target 0.05
    --adaptive-min-frac 1e-8
    --adaptive-grow 4
    --da-phys 5e-6
    --target-crack-extension-um "$TARGET_EXT_UM"
    --crystal-aniso
    --crystal-compete
    --crystal-theta-deg "$CRYSTAL_THETA_DEG"
    --crystal-C11 523e9
    --crystal-C12 203e9
    --crystal-C44 160e9
    --cleave-gamma-aniso 0.3
    --crystal-material w
    --plane-gate-global
    --min-global-forward "$MIN_GLOBAL_FORWARD"
    --max-fronts 1
    --crack-backend adaptive_czm
    --czm-max-angle-error-deg 35
    --j-decomposition cluster
    --mpz-length-um 50
    --mpz-n-bins 80
    --save-snapshots "$SAVE_SNAPSHOTS"
    --snapshot-cols "$SNAPSHOT_COLS"
    --snapshot-by-crack-extension-um "$SNAPSHOT_BY_EXT_UM"
    --out "$OUT"
  )
  if [[ "$GENERATE_SOLVER_PLOTS" == "0" ]]; then
    CMD+=(--no-plots)
  fi

  env \
    CLEAVAGE_HAZARD_MODE=exponential \
    CLEAVAGE_HAZARD_SEED="$SEED" \
    CLEAVAGE_HAZARD_MIN_THRESHOLD=1e-12 \
    CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled \
    CLEAVAGE_EVENT_MIN_FACTOR="$EVENT_MIN_FACTOR" \
    CLEAVAGE_EVENT_MAX_FACTOR="$EVENT_MAX_FACTOR" \
    ARRHENIUS_QUALITY_SUBDIVISION_COUNTS="$QUALITY_SUBDIVISION_COUNTS" \
    ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH="$QUALITY_PARTITION_MAX_DEPTH" \
    ARRHENIUS_QUALITY_PARTITION_MIN_SEGMENT_M="$QUALITY_PARTITION_MIN_SEGMENT_M" \
    "${CMD[@]}" 2>&1 | tee "$LOG"

  echo "[DONE] option=$OPTION T=${T}K seed=$SEED"
}

export -f run_case
export ROOT FAMILY_JSON CAMPAIGN_ROOT TARGET_EXT_UM STEPS DU DT
export SKIP_FINISHED RESTART_INCOMPLETE PRINT_EVERY
export EVENT_MIN_FACTOR EVENT_MAX_FACTOR QUALITY_SUBDIVISION_COUNTS
export QUALITY_PARTITION_MAX_DEPTH QUALITY_PARTITION_MIN_SEGMENT_M
export SNAPSHOT_BY_EXT_UM SAVE_SNAPSHOTS SNAPSHOT_COLS GENERATE_SOLVER_PLOTS
export CRYSTAL_THETA_DEG MIN_GLOBAL_FORWARD

set +e
xargs -P "$MAX_JOBS" -n 4 bash -c 'run_case "$1" "$2" "$3" "$4"' _ < "$CASE_MATRIX"
RUN_STATUS=$?
set -e

python "$ROOT/scripts/analyze_v100517_paper_parameter_campaign.py" \
  "$CAMPAIGN_ROOT" --target-extension-um "$TARGET_EXT_UM" --allow-partial

if [[ "$RUN_STATUS" -ne 0 ]]; then
  echo "ERROR: one or more production cases failed; partial analysis was written." >&2
  exit "$RUN_STATUS"
fi

python "$ROOT/scripts/analyze_v100517_paper_parameter_campaign.py" \
  "$CAMPAIGN_ROOT" --target-extension-um "$TARGET_EXT_UM"

echo "Standalone frozen four-class v10.0.5.19 temperature sweep complete: $CAMPAIGN_ROOT"
