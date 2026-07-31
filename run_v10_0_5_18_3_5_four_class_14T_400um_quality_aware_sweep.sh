#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
RUNNER="$ROOT/run_v10_0_5_18_3_5_four_class_focused_long_growth.sh"

if [[ ! -f "$RUNNER" ]]; then
  echo "ERROR: missing v10.0.5.18.3.5 quality-aware runner: $RUNNER" >&2
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
ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS=${ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS:-8}
ARRHENIUS_QUALITY_AWARE_PARTITIONS=${ARRHENIUS_QUALITY_AWARE_PARTITIONS:-2,4,8}
ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY=${ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY:-0.035}
ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO=${ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO:-0.08}

if [[ -z "${CAMPAIGN_ROOT:-}" ]]; then
  CAMPAIGN_ROOT="$ROOT/runs/v10_0_5_18_3_5_four_class_14T_400um_theta${THETA}_quality_aware_v1"
fi

if pgrep -f 'python.*mode_i_first_passage_v10_0_5_18_3_[12345]' >/dev/null; then
  echo "ERROR: another v10.0.5.18.3.x solver is active" >&2
  ps -axo pid,ppid,etime,%cpu,%mem,state,command |
    grep '[m]ode_i_first_passage_v10_0_5_18_3_' >&2 || true
  exit 1
fi

export ROOT OPTIONS TEMPERATURES THETA MIN_GLOBAL_FORWARD REPLICATE
export TARGET_EXT_UM STEPS MAX_JOBS SKIP_FINISHED BASE_SEED PRINT_EVERY
export GENERATE_SOLVER_PLOTS SAVE_SNAPSHOTS SNAPSHOT_BY_EXT_UM SNAPSHOT_COLS
export ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS ARRHENIUS_QUALITY_AWARE_PARTITIONS
export ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY
export ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO
export CAMPAIGN_ROOT

mkdir -p "$CAMPAIGN_ROOT"
cat > "$CAMPAIGN_ROOT/temperature_sweep_contract.txt" <<EOF
release=10.0.5.18.3.5
runner=run_v10_0_5_18_3_5_four_class_focused_long_growth.sh
sweep_wrapper=run_v10_0_5_18_3_5_four_class_14T_400um_quality_aware_sweep.sh
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
quality_gate_inside_retry_transaction=true
shape_regular_local_tip_patch_bisection=true
quality_aware_patch_levels=$ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS
quality_aware_partition_counts=$ARRHENIUS_QUALITY_AWARE_PARTITIONS
minimum_triangle_quality=$ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY
minimum_immediate_child_area_ratio=$ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO
triangle_quality_floor_relaxed=false
child_area_ratio_floor_relaxed=false
immediate_parent_child_area_ratio_enforced=true
cumulative_child_area_ratio_role=audit_only
exact_stochastic_event_endpoint_preserved=true
exact_selected_crack_direction_preserved=true
path_aware_growth_envelope=true
EOF

bash "$RUNNER"
