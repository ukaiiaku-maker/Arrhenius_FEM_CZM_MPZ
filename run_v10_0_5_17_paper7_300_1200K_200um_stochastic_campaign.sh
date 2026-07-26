#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PARAMETER_SOURCE_ROOT=${PARAMETER_SOURCE_ROOT:-$ROOT/runtime_inputs/v10_0_5_17_frozen_pf_inputs}
FAMILY_JSON=${FAMILY_JSON:-$PARAMETER_SOURCE_ROOT/signed_kernel/v10_2_14_active_only_campaign_family.json}
FAMILY_JSON_SHA256_REQUIRED=${FAMILY_JSON_SHA256_REQUIRED:-a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_17_paper7_300_1200K_200um_stochastic_pf_parity_v1}

TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 850 900 950 1000 1050 1100 1150 1200"}
TARGET_EXT_UM=${TARGET_EXT_UM:-200}
STEPS=${STEPS:-100000}
DU=${DU:-2e-5}
DT=${DT:-840}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
PRINT_EVERY=${PRINT_EVERY:-100}
BASE_HAZARD_SEED=${BASE_HAZARD_SEED:-1720}
EVENT_MIN_FACTOR=${EVENT_MIN_FACTOR:-0.5}
EVENT_MAX_FACTOR=${EVENT_MAX_FACTOR:-4.0}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-10}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-20}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-1}
INCLUDE_CONTROL=${INCLUDE_CONTROL:-0}

PRIMARY_OPTIONS=(
  v913_paper_peak01_0242980_persistent_sites
  v913_paper_peak02_0127508_persistent_sites
  v913_paper_peak03_0115460_persistent_sites
  v913_paper_dbtt01_0202500_persistent_sites
  v913_paper_dbtt02_0088403_persistent_sites
  v913_paper_weakT01_0257068_persistent_sites
  v913_paper_ceramic01_0189364_persistent_sites
)
CONTROL_OPTION=v913_paper_control01_0086420_persistent_sites

MANIFEST="$PARAMETER_SOURCE_ROOT/frozen_input_manifest_v10_0_5_17.json"
CATALOG="$PARAMETER_SOURCE_ROOT/frozen_parameter_catalog_v10_0_5_17.json"
REQUIRED_INPUTS=(
  "$MANIFEST"
  "$CATALOG"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/data/materials/v10_2_25_v913_paper_campaign_registry.csv"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/data/materials/v10_2_25_v913_paper_campaign_selection.json"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/sharp_front_v10_2_25.py"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/sharp_front_v10_2_25_audited.py"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/data/materials/v10_2_26_v913_weakT_ceramic_registry.csv"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/data/materials/v10_2_26_v913_weakT_ceramic_selection.json"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/sharp_front_v10_2_26.py"
  "$PARAMETER_SOURCE_ROOT/arrhenius_fracture/sharp_front_v10_2_26_audited.py"
  "$FAMILY_JSON"
)
for path in "${REQUIRED_INPUTS[@]}"; do
  if [[ ! -f "$path" ]]; then
    echo "ERROR: missing local frozen FEM/CZM input: $path" >&2
    exit 1
  fi
done

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

python - "$MANIFEST" "$FAMILY_JSON_SHA256" <<'PY'
import json
import pathlib
import sys
manifest = json.loads(pathlib.Path(sys.argv[1]).read_text())
observed = sys.argv[2]
errors = []
if manifest.get("solver_or_constitutive_files_modified") is not False:
    errors.append("vendor manifest reports solver or constitutive changes")
if manifest.get("parameter_values_reconstructed") is not False:
    errors.append("vendor manifest reports reconstructed parameter values")
if manifest.get("kernel_family_vendored_byte_for_byte") is not True:
    errors.append("kernel family was not vendored byte-for-byte")
if manifest.get("kernel_family_sha256") != observed:
    errors.append("vendor manifest kernel SHA does not match local family")
if errors:
    raise SystemExit("; ".join(errors))
PY

OPTIONS=("${PRIMARY_OPTIONS[@]}")
if [[ "$INCLUDE_CONTROL" == "1" ]]; then
  OPTIONS+=("$CONTROL_OPTION")
elif [[ "$INCLUDE_CONTROL" != "0" ]]; then
  echo "ERROR: INCLUDE_CONTROL must be 0 or 1" >&2
  exit 1
fi

mkdir -p "$CAMPAIGN_ROOT"
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
entry = {
    "v913_paper_peak01_0242980_persistent_sites": "v10.2.25",
    "v913_paper_peak02_0127508_persistent_sites": "v10.2.25",
    "v913_paper_peak03_0115460_persistent_sites": "v10.2.25",
    "v913_paper_dbtt01_0202500_persistent_sites": "v10.2.25",
    "v913_paper_dbtt02_0088403_persistent_sites": "v10.2.25",
    "v913_paper_control01_0086420_persistent_sites": "v10.2.25",
    "v913_paper_weakT01_0257068_persistent_sites": "v10.2.26",
    "v913_paper_ceramic01_0189364_persistent_sites": "v10.2.26",
}
rows = []
seen = set()
for option in options:
    if option not in entry:
        raise SystemExit(f"unknown paper option: {option}")
    for temperature in temperatures:
        token = f"v10.0.5.17|{base}|{option}|{temperature}".encode()
        seed = int.from_bytes(hashlib.sha256(token).digest()[:4], "big")
        if seed in seen:
            raise SystemExit(f"seed collision for {option} at {temperature} K")
        seen.add(seed)
        rows.append((entry[option], option, temperature, seed))
Path(os.environ["CASE_MATRIX"]).write_text(
    "".join(f"{a}\t{b}\t{c}\t{d}\n" for a, b, c, d in rows)
)
with Path(os.environ["SEED_MAP"]).open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["parameter_entry", "parameter_option", "temperature_K", "hazard_seed"])
    writer.writerows(rows)
PY

cat > "$CAMPAIGN_ROOT/campaign_configuration.txt" <<EOF
release=10.0.5.17
entry=arrhenius_fracture.mode_i_first_passage_v10_0_5_17_paper_parameter_campaign
parameter_source_root=$PARAMETER_SOURCE_ROOT
external_PF_runtime_dependency=false
parameter_source_vendored_into_FEM_CZM_checkout=true
parameter_options=${OPTIONS[*]}
include_rehardening_control=$INCLUDE_CONTROL
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
kernel_family=$FAMILY_JSON
kernel_family_sha256=$FAMILY_JSON_SHA256
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

  PLOT_ARGS=()
  if [[ "$GENERATE_SOLVER_PLOTS" == "0" ]]; then
    PLOT_ARGS+=(--no-plots)
  fi

  env \
    CLEAVAGE_HAZARD_MODE=exponential \
    CLEAVAGE_HAZARD_SEED="$SEED" \
    CLEAVAGE_HAZARD_MIN_THRESHOLD=1e-12 \
    CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled \
    CLEAVAGE_EVENT_MIN_FACTOR="$EVENT_MIN_FACTOR" \
    CLEAVAGE_EVENT_MAX_FACTOR="$EVENT_MAX_FACTOR" \
    python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_17_paper_parameter_campaign \
      --parameter-source-root "$PARAMETER_SOURCE_ROOT" \
      --parameter-entry "$ENTRY" \
      --parameter-option "$OPTION" \
      --signed-kernel-family "$FAMILY_JSON" \
      --tip-refinement-radius-um 330 \
      --selected-cluster-J-outer-um 240 \
      --local-J-outer-um 100 \
      --mode 2d \
      --bulk-plasticity-mode tip_only \
      --temperatures "$T" \
      --steps "$STEPS" \
      --nx 36 --ny 72 \
      --tip-h-fine 2.5e-6 --tip-ratio 1.15 \
      --dU "$DU" --dt "$DT" \
      --n-stagger 1 \
      --print-every "$PRINT_EVERY" \
      --adaptive-events \
      --adaptive-event-target 0.05 \
      --adaptive-min-frac 1e-8 \
      --adaptive-grow 4 \
      --da-phys 5e-6 \
      --target-crack-extension-um "$TARGET_EXT_UM" \
      --crystal-aniso \
      --crystal-compete \
      --crystal-theta-deg 45 \
      --crystal-C11 523e9 \
      --crystal-C12 203e9 \
      --crystal-C44 160e9 \
      --cleave-gamma-aniso 0.3 \
      --crystal-material w \
      --max-fronts 1 \
      --crack-backend adaptive_czm \
      --czm-max-angle-error-deg 35 \
      --j-decomposition cluster \
      --mpz-length-um 50 \
      --mpz-n-bins 80 \
      --save-snapshots "$SAVE_SNAPSHOTS" \
      --snapshot-cols "$SNAPSHOT_COLS" \
      --snapshot-by-crack-extension-um "$SNAPSHOT_BY_EXT_UM" \
      "${PLOT_ARGS[@]}" \
      --out "$OUT" \
      2>&1 | tee "$LOG"

  echo "[DONE] option=$OPTION T=${T}K seed=$SEED"
}

export -f run_case
export ROOT PARAMETER_SOURCE_ROOT FAMILY_JSON CAMPAIGN_ROOT TARGET_EXT_UM STEPS DU DT
export SKIP_FINISHED PRINT_EVERY EVENT_MIN_FACTOR EVENT_MAX_FACTOR
export SNAPSHOT_BY_EXT_UM SAVE_SNAPSHOTS SNAPSHOT_COLS GENERATE_SOLVER_PLOTS

xargs -P "$MAX_JOBS" -n 4 bash -c 'run_case "$1" "$2" "$3" "$4"' _ < "$CASE_MATRIX"

python "$ROOT/scripts/analyze_v100517_paper_parameter_campaign.py" \
  "$CAMPAIGN_ROOT" \
  --target-extension-um "$TARGET_EXT_UM"

echo "Campaign and final analysis complete: $CAMPAIGN_ROOT"
