#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
RUNNER="$ROOT/run_v10_0_5_18_3_4_four_class_focused_long_growth.sh"

if [[ ! -f "$RUNNER" ]]; then
  echo "ERROR: missing v10.0.5.18.3.4 path-envelope runner: $RUNNER" >&2
  exit 1
fi

OPTIONS=${OPTIONS:-"v913_paper_peak01_0242980_persistent_sites v913_paper_dbtt01_0202500_persistent_sites v913_paper_weakT01_0129902_persistent_sites v913_paper_ceramic01_0077080_persistent_sites"}
TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 850 900 950 1000 1050 1100 1150 1200"}
THETA=${THETA:-30}
MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}
REPLICATE=${REPLICATE:-0}
TARGET_EXT_UM=${TARGET_EXT_UM:-400}
STEPS=${STEPS:-1000000}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
BASE_SEED=${BASE_SEED:-1720}
PRINT_EVERY=${PRINT_EVERY:-50}
GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-0}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-8}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-4}
ARRHENIUS_CORRIDOR_GUARD_UM=${ARRHENIUS_CORRIDOR_GUARD_UM:-10}
ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ=${ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ:-0.25}
ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ=${ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ:-0.245}
ARRHENIUS_ENVELOPE_BUFFER_LPZ=${ARRHENIUS_ENVELOPE_BUFFER_LPZ:-1.0}
ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES=${ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES:-181}
ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES=${ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES:-4}
ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR=${ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR:-0.9}
ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ=${ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ:-0.5}
ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS=${ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS:-16}
EMISSION_EVENT_HORIZON_FACTOR=${EMISSION_EVENT_HORIZON_FACTOR:-1.25}
EMISSION_MAX_EVENTS_PER_HALF_STEP=${EMISSION_MAX_EVENTS_PER_HALF_STEP:-10000}

if [[ -z "${CAMPAIGN_ROOT:-}" ]]; then
  CAMPAIGN_ROOT="$ROOT/runs/v10_0_5_18_3_4_four_class_14T_400um_theta${THETA}_path_envelope_v1"
fi

if pgrep -f 'python.*mode_i_first_passage_v10_0_5_18_3_[1234]' >/dev/null; then
  echo "ERROR: another v10.0.5.18.3.x solver is active" >&2
  ps -axo pid,ppid,etime,%cpu,%mem,state,command |
    grep '[m]ode_i_first_passage_v10_0_5_18_3_' >&2 || true
  exit 1
fi

export ROOT OPTIONS TEMPERATURES THETA MIN_GLOBAL_FORWARD REPLICATE
export TARGET_EXT_UM STEPS MAX_JOBS SKIP_FINISHED BASE_SEED PRINT_EVERY
export GENERATE_SOLVER_PLOTS SAVE_SNAPSHOTS SNAPSHOT_BY_EXT_UM SNAPSHOT_COLS
export ARRHENIUS_CORRIDOR_GUARD_UM ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ
export ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ ARRHENIUS_ENVELOPE_BUFFER_LPZ
export ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES
export ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR
export ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS
export EMISSION_EVENT_HORIZON_FACTOR EMISSION_MAX_EVENTS_PER_HALF_STEP
export CAMPAIGN_ROOT

mkdir -p "$CAMPAIGN_ROOT"
cat > "$CAMPAIGN_ROOT/temperature_sweep_contract.txt" <<EOF
release=10.0.5.18.3.4
runner=run_v10_0_5_18_3_4_four_class_focused_long_growth.sh
sweep_wrapper=run_v10_0_5_18_3_4_four_class_14T_400um_path_envelope_sweep.sh
parameter_options=$OPTIONS
temperatures_K=$TEMPERATURES
theta_deg=$THETA
minimum_global_forward_cosine=$MIN_GLOBAL_FORWARD
replicate=$REPLICATE
target_extension_um=$TARGET_EXT_UM
steps_ceiling=$STEPS
max_jobs=$MAX_JOBS
skip_finished=$SKIP_FINISHED
base_seed=$BASE_SEED
case_seed_namespace=v10.0.5.18.3
save_snapshots=$SAVE_SNAPSHOTS
snapshot_by_extension_um=$SNAPSHOT_BY_EXT_UM
snapshot_columns=$SNAPSHOT_COLS
corridor_guard_um=$ARRHENIUS_CORRIDOR_GUARD_UM
maximum_corridor_h_tip_over_L_pz=$ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ
envelope_support_h_over_L_pz=$ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ
envelope_buffer_L_pz=$ARRHENIUS_ENVELOPE_BUFFER_LPZ
envelope_direction_samples=$ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES
emission_event_horizon_factor=$EMISSION_EVENT_HORIZON_FACTOR
emission_max_events_per_half_step=$EMISSION_MAX_EVENTS_PER_HALF_STEP
single_trial_event_horizon_active=true
state_dependent_log_rate_rejection_active=false
PF_transport_operator_preserved=true
emission_events_batched=false
path_aware_growth_envelope=true
directional_support_basis=complete_global_forward_gate_conservative
full_reachable_envelope_process_zone_resolved=true
production_physical_refinement_provider_required=true
initial_physical_refinement_radius_um=330
swept_physical_refinement_active=false
child_area_ratio_floor_relaxed=false
EOF

bash "$RUNNER"
