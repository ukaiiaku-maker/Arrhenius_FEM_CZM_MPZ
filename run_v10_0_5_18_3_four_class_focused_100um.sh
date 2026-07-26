#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PARAMETER_SOURCE_ROOT=${PARAMETER_SOURCE_ROOT:-$ROOT/runtime_inputs/v10_0_5_18_four_class}
FAMILY_JSON=${FAMILY_JSON:-$ROOT/runtime_inputs/v10_0_5_17_frozen_pf_inputs/signed_kernel/v10_2_14_active_only_campaign_family.json}
FAMILY_JSON_SHA256_REQUIRED=${FAMILY_JSON_SHA256_REQUIRED:-a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_18_3_weakT_ceramic_5T_100um_load_ramp_theta${THETA:-30}_v1}

OPTIONS=${OPTIONS:-"v913_paper_weakT01_0129902_persistent_sites v913_paper_ceramic01_0077080_persistent_sites"}
TEMPERATURES=${TEMPERATURES:-"300 600 900 1100 1200"}
THETA=${THETA:-30}
REPLICATE=${REPLICATE:-0}
TARGET_EXT_UM=${TARGET_EXT_UM:-100}
STEPS=${STEPS:-100000}
DU=${DU:-2e-5}
DT=${DT:-840}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
PRINT_EVERY=${PRINT_EVERY:-100}
BASE_SEED=${BASE_SEED:-1720}
EVENT_MIN_FACTOR=${EVENT_MIN_FACTOR:-0.5}
EVENT_MAX_FACTOR=${EVENT_MAX_FACTOR:-4.0}
ADAPTIVE_EVENT_TARGET=${ADAPTIVE_EVENT_TARGET:-0.05}
EMISSION_MAX_ACTION_SUBSTEP=${EMISSION_MAX_ACTION_SUBSTEP:-0.05}
EMISSION_MAX_EVENTS_PER_HALF_STEP=${EMISSION_MAX_EVENTS_PER_HALF_STEP:-10000}
EMISSION_RAMP_MAX_LOG_RATE_CHANGE=${EMISSION_RAMP_MAX_LOG_RATE_CHANGE:-0.5}
EMISSION_RAMP_ACTION_FLOOR=${EMISSION_RAMP_ACTION_FLOOR:-1e-6}
TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS=${TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS:-8}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-10}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-10}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-1}
PYTHONUNBUFFERED=${PYTHONUNBUFFERED:-1}

REGISTRY="$PARAMETER_SOURCE_ROOT/v10_2_27_v913_four_class_paper_registry.csv"
TRANSFER_MANIFEST="$PARAMETER_SOURCE_ROOT/four_class_transfer_manifest_v10_0_5_18.json"
for path in "$REGISTRY" "$TRANSFER_MANIFEST" "$FAMILY_JSON"; do
  if [[ ! -f "$path" ]]; then
    echo "ERROR: required local FEM/CZM input is missing: $path" >&2
    exit 1
  fi
done

FAMILY_JSON=$(cd "$(dirname "$FAMILY_JSON")" && pwd)/$(basename "$FAMILY_JSON")
FAMILY_JSON_SHA256=$(python - "$FAMILY_JSON" <<'PY'
import hashlib
import pathlib
import sys
print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)
if [[ "$FAMILY_JSON_SHA256" != "$FAMILY_JSON_SHA256_REQUIRED" ]]; then
  echo "ERROR: signed-kernel family SHA-256 mismatch" >&2
  echo "Expected: $FAMILY_JSON_SHA256_REQUIRED" >&2
  echo "Observed: $FAMILY_JSON_SHA256" >&2
  exit 1
fi

mkdir -p "$CAMPAIGN_ROOT"
CASE_MATRIX="$CAMPAIGN_ROOT/case_matrix.tsv"
SEED_MAP="$CAMPAIGN_ROOT/stochastic_seed_map.csv"

OPTIONS_TEXT=$(printf '%s\n' $OPTIONS)
OPTIONS_TEXT="$OPTIONS_TEXT" \
TEMPERATURES="$TEMPERATURES" \
BASE_SEED="$BASE_SEED" \
REPLICATE="$REPLICATE" \
THETA="$THETA" \
CASE_MATRIX="$CASE_MATRIX" \
SEED_MAP="$SEED_MAP" \
python - <<'PY'
import csv
import hashlib
import os
from pathlib import Path

options = [line.strip() for line in os.environ["OPTIONS_TEXT"].splitlines() if line.strip()]
temperatures = [int(value) for value in os.environ["TEMPERATURES"].split()]
base = int(os.environ["BASE_SEED"])
replicate = int(os.environ["REPLICATE"])
theta = float(os.environ["THETA"])
allowed = {
    "v913_paper_peak01_0242980_persistent_sites",
    "v913_paper_dbtt01_0202500_persistent_sites",
    "v913_paper_weakT01_0129902_persistent_sites",
    "v913_paper_ceramic01_0077080_persistent_sites",
}

def seed_for(option: str, temperature: int) -> int:
    token = (
        f"v10.0.5.18.3|base={base}|option={option}|T={temperature}|"
        f"replicate={replicate}"
    ).encode("utf-8")
    return int.from_bytes(hashlib.sha256(token).digest()[:4], "big")

def stream_seed(case_seed: int, channel: int) -> int:
    token = (
        f"v10.0.5.18|case_seed={case_seed}|stream=signed_emission|"
        f"channel={channel}"
    ).encode("utf-8")
    return int.from_bytes(hashlib.sha256(token).digest()[:8], "big")

rows = []
seen = set()
for option in options:
    if option not in allowed:
        raise SystemExit(f"unknown four-class option: {option}")
    for temperature in temperatures:
        seed = seed_for(option, temperature)
        if seed in seen:
            raise SystemExit(f"case-seed collision for {option} at {temperature} K")
        seen.add(seed)
        rows.append(
            (
                option,
                temperature,
                seed,
                stream_seed(seed, 0),
                stream_seed(seed, 1),
                replicate,
                theta,
            )
        )

Path(os.environ["CASE_MATRIX"]).write_text(
    "".join(f"{option}\t{temperature}\t{seed}\n" for option, temperature, seed, *_ in rows)
)
with Path(os.environ["SEED_MAP"]).open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(
        [
            "parameter_option",
            "temperature_K",
            "cleavage_case_seed",
            "emission_channel_0_seed",
            "emission_channel_1_seed",
            "replicate",
            "theta_deg",
        ]
    )
    writer.writerows(rows)
PY

cat > "$CAMPAIGN_ROOT/campaign_configuration.txt" <<EOF
release=10.0.5.18.3
entry=arrhenius_fracture.mode_i_first_passage_v10_0_5_18_3_four_class_load_ramp_stochastic_emission
parameter_source_root=$PARAMETER_SOURCE_ROOT
parameter_options=$OPTIONS
temperatures_K=$TEMPERATURES
theta_deg=$THETA
theta_is_runtime_input=true
replicate=$REPLICATE
target_extension_um=$TARGET_EXT_UM
base_seed=$BASE_SEED
cleavage_hazard=exponential_integrated_action
cleavage_event_length=threshold_correlated_mean_preserving
emission_hazard=independent_exponential_integrated_action_per_signed_channel
emission_event_packet=one_activation_times_mechanical_line_conversion
emission_accepted_update=exact_event_localized
emission_post_onset_mean_field_burst=false
global_FEM_deterministic_emission_state_predictor=false
local_linear_K_ramp_kinetics=true
emission_scalar_directional_factor_applied_to_opening_K=false
emission_signed_channel_factors_applied_once=true
adaptive_event_target=$ADAPTIVE_EVENT_TARGET
emission_max_action_substep=$EMISSION_MAX_ACTION_SUBSTEP
emission_ramp_max_log_rate_change=$EMISSION_RAMP_MAX_LOG_RATE_CHANGE
emission_ramp_action_floor=$EMISSION_RAMP_ACTION_FLOOR
two_channel_probe_fallback_max_calls=$TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS
adaptive_CZM_atomic_partition_retry=true
continuous_fractional_MPZ_translation=true
microstructure_advance_precedes_cohesive_checkpoint=true
kernel_family=$FAMILY_JSON
kernel_family_sha256=$FAMILY_JSON_SHA256
max_jobs=$MAX_JOBS
EOF

run_case() {
  set -euo pipefail
  local OPTION="$1"
  local T="$2"
  local CASE_SEED="$3"
  local TAG
  printf -v TAG "%04d" "$T"
  local OUT="$CAMPAIGN_ROOT/$OPTION/T${TAG}K"
  local MANIFEST="$OUT/persistent_site_production_manifest_v10_0_5_18_3.json"
  local LOG="$OUT/console.log"

  if [[ "$SKIP_FINISHED" == "1" && -f "$MANIFEST" ]] &&
     grep -q '"run_completed_without_exception": true' "$MANIFEST"; then
    echo "[SKIP] option=$OPTION T=${T}K"
    return 0
  fi

  rm -rf "$OUT"
  mkdir -p "$OUT"
  echo "[START] option=$OPTION T=${T}K theta=${THETA} seed=$CASE_SEED"

  PLOT_ARGS=()
  if [[ "$GENERATE_SOLVER_PLOTS" == "0" ]]; then
    PLOT_ARGS+=(--no-plots)
  fi

  env \
    PYTHONUNBUFFERED="$PYTHONUNBUFFERED" \
    CLEAVAGE_HAZARD_MODE=exponential \
    CLEAVAGE_HAZARD_SEED="$CASE_SEED" \
    CLEAVAGE_HAZARD_MIN_THRESHOLD=1e-12 \
    CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled \
    CLEAVAGE_EVENT_MIN_FACTOR="$EVENT_MIN_FACTOR" \
    CLEAVAGE_EVENT_MAX_FACTOR="$EVENT_MAX_FACTOR" \
    EMISSION_HAZARD_SEED="$CASE_SEED" \
    EMISSION_HAZARD_MIN_THRESHOLD=1e-12 \
    EMISSION_MAX_ACTION_SUBSTEP="$EMISSION_MAX_ACTION_SUBSTEP" \
    EMISSION_MAX_EVENTS_PER_HALF_STEP="$EMISSION_MAX_EVENTS_PER_HALF_STEP" \
    EMISSION_RAMP_MAX_LOG_RATE_CHANGE="$EMISSION_RAMP_MAX_LOG_RATE_CHANGE" \
    EMISSION_RAMP_ACTION_FLOOR="$EMISSION_RAMP_ACTION_FLOOR" \
    TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS="$TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS" \
    python -u -m arrhenius_fracture.mode_i_first_passage_v10_0_5_18_3_four_class_load_ramp_stochastic_emission \
      --parameter-source-root "$PARAMETER_SOURCE_ROOT" \
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
      --adaptive-event-target "$ADAPTIVE_EVENT_TARGET" \
      --adaptive-min-frac 1e-8 \
      --adaptive-grow 4 \
      --da-phys 5e-6 \
      --target-crack-extension-um "$TARGET_EXT_UM" \
      --crystal-aniso \
      --crystal-compete \
      --crystal-theta-deg "$THETA" \
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

  echo "[DONE] option=$OPTION T=${T}K theta=${THETA} seed=$CASE_SEED"
}

export -f run_case
export ROOT PARAMETER_SOURCE_ROOT FAMILY_JSON CAMPAIGN_ROOT TARGET_EXT_UM STEPS DU DT
export SKIP_FINISHED PRINT_EVERY EVENT_MIN_FACTOR EVENT_MAX_FACTOR THETA
export ADAPTIVE_EVENT_TARGET EMISSION_MAX_ACTION_SUBSTEP EMISSION_MAX_EVENTS_PER_HALF_STEP
export EMISSION_RAMP_MAX_LOG_RATE_CHANGE EMISSION_RAMP_ACTION_FLOOR
export TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS PYTHONUNBUFFERED
export SNAPSHOT_BY_EXT_UM SAVE_SNAPSHOTS SNAPSHOT_COLS GENERATE_SOLVER_PLOTS

xargs -P "$MAX_JOBS" -n 3 bash -c 'run_case "$1" "$2" "$3"' _ < "$CASE_MATRIX"

echo "Focused v10.0.5.18.3 validation complete: $CAMPAIGN_ROOT"
