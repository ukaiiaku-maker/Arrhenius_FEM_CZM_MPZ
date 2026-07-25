#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PFROOT=${PFROOT:-/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1}
PF_BRANCH_REQUIRED=${PF_BRANCH_REQUIRED:-v10.2.22-physical-front-width-top5-dbtt-screen}
PF_COMMIT_PREFIX_REQUIRED=${PF_COMMIT_PREFIX_REQUIRED:-0a340f6}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_17_v10227_paper4_anisotropic_300_1200K_1000um_stochastic_pf_parity_v1}

ATLAS_REL=runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json
FAMILY_JSON=${FAMILY_JSON:-}
FAMILY_JSON_SOURCE=explicit_environment
if [[ -z "$FAMILY_JSON" ]]; then
  FAMILY_JSON_SOURCE=unresolved
  ATLAS_CANDIDATES=(
    "$PFROOT/$ATLAS_REL"
    "$(dirname "$ROOT")/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/$ATLAS_REL"
  )
  for candidate in "${ATLAS_CANDIDATES[@]}"; do
    if [[ -f "$candidate" ]]; then
      FAMILY_JSON="$candidate"
      FAMILY_JSON_SOURCE=validated_runtime_artifact
      break
    fi
  done
fi
if [[ -z "$FAMILY_JSON" || ! -f "$FAMILY_JSON" ]]; then
  echo "ERROR: missing audited v10.2.14 active-only signed-kernel family." >&2
  echo "Provide FAMILY_JSON=/absolute/path/to/$ATLAS_REL" >&2
  exit 1
fi
FAMILY_JSON=$(cd "$(dirname "$FAMILY_JSON")" && pwd)/$(basename "$FAMILY_JSON")

TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 850 900 950 1000 1050 1100 1150 1200"}
TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
STEPS=${STEPS:-1000000}
DU=${DU:-2e-5}
DT=${DT:-840}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
PRINT_EVERY=${PRINT_EVERY:-100}
BASE_HAZARD_SEED=${BASE_HAZARD_SEED:-1720}
EVENT_MIN_FACTOR=${EVENT_MIN_FACTOR:-0.5}
EVENT_MAX_FACTOR=${EVENT_MAX_FACTOR:-4.0}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-20}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-1}
CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-30}
MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}
RUN_J_SCALE_DIAGNOSTIC=${RUN_J_SCALE_DIAGNOSTIC:-1}
J_DIAGNOSTIC_STRICT=${J_DIAGNOSTIC_STRICT:-0}

for VALUE_NAME in GENERATE_SOLVER_PLOTS RUN_J_SCALE_DIAGNOSTIC J_DIAGNOSTIC_STRICT; do
  VALUE=${!VALUE_NAME}
  case "$VALUE" in
    0|1) ;;
    *) echo "ERROR: $VALUE_NAME must be 0 or 1" >&2; exit 1 ;;
  esac
done

OPTIONS=(
  v913_paper_peak01_0242980_persistent_sites
  v913_paper_dbtt01_0202500_persistent_sites
  v913_paper_weakT01_0129902_persistent_sites
  v913_paper_ceramic01_0077080_persistent_sites
)

if [[ ! -d "$PFROOT/.git" ]]; then
  echo "ERROR: PFROOT is not a Git clone: $PFROOT" >&2
  exit 1
fi
PF_BRANCH=$(git -C "$PFROOT" branch --show-current)
PF_COMMIT=$(git -C "$PFROOT" rev-parse HEAD)
if [[ "$PF_BRANCH" != "$PF_BRANCH_REQUIRED" ]]; then
  echo "ERROR: PFROOT is on '$PF_BRANCH'; required '$PF_BRANCH_REQUIRED'" >&2
  exit 1
fi
if [[ "$PF_COMMIT" != "$PF_COMMIT_PREFIX_REQUIRED"* ]]; then
  echo "ERROR: PFROOT commit '$PF_COMMIT' does not match canonical prefix '$PF_COMMIT_PREFIX_REQUIRED'" >&2
  exit 1
fi

REQUIRED_PF_FILES=(
  "$PFROOT/arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"
  "$PFROOT/arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_selection.json"
  "$PFROOT/arrhenius_fracture/sharp_front_v10_2_27.py"
  "$PFROOT/arrhenius_fracture/sharp_front_v10_2_27_audited.py"
  "$FAMILY_JSON"
)
for path in "${REQUIRED_PF_FILES[@]}"; do
  if [[ ! -f "$path" ]]; then
    echo "ERROR: missing audited PF input: $path" >&2
    exit 1
  fi
done

mkdir -p "$CAMPAIGN_ROOT"
python "$ROOT/scripts/validate_v100517_v10227_four_class_inputs.py" \
  --pf-repo-root "$PFROOT" \
  --out "$CAMPAIGN_ROOT/audited_parameter_catalog.json" \
  "${OPTIONS[@]}"

FAMILY_JSON_SHA256=$(python - "$FAMILY_JSON" <<'PY'
from pathlib import Path
import hashlib
import sys
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)

echo "Using PF v10.2.27 four-class registry at commit $PF_COMMIT"
echo "Using signed-kernel family: $FAMILY_JSON"
echo "  source: $FAMILY_JSON_SOURCE"
echo "  sha256: $FAMILY_JSON_SHA256"
echo "Parameter set: paper4-v10.2.27 (${#OPTIONS[@]} options)"
echo "Crystallographic path selection: theta=${CRYSTAL_THETA_DEG} deg, global forward gate=${MIN_GLOBAL_FORWARD}"

if [[ "$RUN_J_SCALE_DIAGNOSTIC" == "1" ]]; then
  JDIAG_ARGS=(
    python "$ROOT/scripts/diagnose_v100517_j_scale.py"
    --out "$CAMPAIGN_ROOT/preflight_J_scale_diagnostic"
    --theta-deg "$CRYSTAL_THETA_DEG"
  )
  if [[ "$J_DIAGNOSTIC_STRICT" == "1" ]]; then
    JDIAG_ARGS+=(--strict)
  fi
  "${JDIAG_ARGS[@]}"
fi

CASE_MATRIX="$CAMPAIGN_ROOT/case_matrix.tsv"
SEED_MAP="$CAMPAIGN_ROOT/parameter_temperature_seed_map.csv"
OPTIONS_TEXT=$(printf '%s\n' "${OPTIONS[@]}")
OPTIONS_TEXT="$OPTIONS_TEXT" TEMPERATURES="$TEMPERATURES" BASE_HAZARD_SEED="$BASE_HAZARD_SEED" \
CASE_MATRIX="$CASE_MATRIX" SEED_MAP="$SEED_MAP" python - <<'PY'
import csv
import hashlib
import os
from pathlib import Path

options = [line.strip() for line in os.environ["OPTIONS_TEXT"].splitlines() if line.strip()]
temperatures = [int(value) for value in os.environ["TEMPERATURES"].split()]
base = int(os.environ["BASE_HAZARD_SEED"])
rows = []
seen = set()
for option in options:
    for temperature in temperatures:
        token = f"v10.0.5.17|{base}|{option}|{temperature}".encode()
        seed = int.from_bytes(hashlib.sha256(token).digest()[:4], "big")
        if seed in seen:
            raise SystemExit(f"seed collision for {option} at {temperature} K")
        seen.add(seed)
        rows.append(("v10.2.27", option, temperature, seed))
Path(os.environ["CASE_MATRIX"]).write_text(
    "".join(f"{a}\t{b}\t{c}\t{d}\n" for a, b, c, d in rows)
)
with Path(os.environ["SEED_MAP"]).open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["parameter_entry", "parameter_option", "temperature_K", "hazard_seed"])
    writer.writerows(rows)
PY

cat > "$CAMPAIGN_ROOT/campaign_configuration.txt" <<EOF
release=10.0.5.17-v10.2.27-parameters
entry=arrhenius_fracture.mode_i_first_passage_v10_0_5_17_v10227_four_class
PF_repo_root=$PFROOT
PF_branch_required=$PF_BRANCH_REQUIRED
PF_commit_prefix_required=$PF_COMMIT_PREFIX_REQUIRED
PF_commit=$PF_COMMIT
parameter_set=paper4-v10.2.27
parameter_options=${OPTIONS[*]}
temperatures_K=$TEMPERATURES
target_extension_um=$TARGET_EXT_UM
base_hazard_seed=$BASE_HAZARD_SEED
hazard_seed_formula=sha256(v10.0.5.17|base|parameter_option|temperature_K)_first_32_bits
cleavage_hazard_mode=minus_log_uniform_exponential_integrated_action
cleavage_event_length_mode=threshold_scaled_same_Xi_mean_preserving
event_min_factor=$EVENT_MIN_FACTOR
event_max_factor=$EVENT_MAX_FACTOR
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
kernel_family_source=$FAMILY_JSON_SOURCE
kernel_family_sha256=$FAMILY_JSON_SHA256
run_J_scale_diagnostic=$RUN_J_SCALE_DIAGNOSTIC
J_diagnostic_strict=$J_DIAGNOSTIC_STRICT
max_jobs=$MAX_JOBS
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

  rm -rf "$OUT"
  mkdir -p "$OUT"
  echo "[START] option=$OPTION entry=$ENTRY T=${T}K seed=$SEED"

  local -a CMD
  CMD=(
    python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_17_v10227_four_class
    --pf-repo-root "$PFROOT"
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
    "${CMD[@]}" 2>&1 | tee "$LOG"

  python "$ROOT/scripts/analyze_v100517_paper_parameter_campaign.py" \
    "$CAMPAIGN_ROOT" --target-extension-um "$TARGET_EXT_UM" --allow-partial \
    --lock-file "$CAMPAIGN_ROOT/.analysis.lock" || true

  echo "[DONE] option=$OPTION T=${T}K seed=$SEED"
}

export -f run_case
export ROOT PFROOT FAMILY_JSON CAMPAIGN_ROOT TARGET_EXT_UM STEPS DU DT
export SKIP_FINISHED PRINT_EVERY EVENT_MIN_FACTOR EVENT_MAX_FACTOR
export SNAPSHOT_BY_EXT_UM SAVE_SNAPSHOTS SNAPSHOT_COLS GENERATE_SOLVER_PLOTS
export CRYSTAL_THETA_DEG MIN_GLOBAL_FORWARD

set +e
xargs -P "$MAX_JOBS" -n 4 bash -c 'run_case "$1" "$2" "$3" "$4"' _ < "$CASE_MATRIX"
RUN_STATUS=$?
set -e

python "$ROOT/scripts/analyze_v100517_paper_parameter_campaign.py" \
  "$CAMPAIGN_ROOT" --target-extension-um "$TARGET_EXT_UM" --allow-partial

if [[ "$RUN_STATUS" -ne 0 ]]; then
  echo "ERROR: one or more cases failed; partial analysis was written." >&2
  exit "$RUN_STATUS"
fi

python "$ROOT/scripts/analyze_v100517_paper_parameter_campaign.py" \
  "$CAMPAIGN_ROOT" --target-extension-um "$TARGET_EXT_UM"

echo "Campaign and final analysis complete: $CAMPAIGN_ROOT"
