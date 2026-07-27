#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_2_four_class_focused_100um.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.2 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
export ARRHENIUS_CORRIDOR_GUARD_UM=${ARRHENIUS_CORRIDOR_GUARD_UM:-10}
export ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM=${ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM:-100}
export ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ=${ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ:-0.25}
export EMISSION_EVENT_HORIZON_FACTOR=${EMISSION_EVENT_HORIZON_FACTOR:-1.25}
export EMISSION_MAX_EVENTS_PER_HALF_STEP=${EMISSION_MAX_EVENTS_PER_HALF_STEP:-10000}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051833_long_growth_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission/mode_i_first_passage_v10_0_5_18_3_3_four_class_long_growth_corridor/g' \
  -e 's/persistent_site_production_manifest_v10_0_5_18_3_2.json/persistent_site_production_manifest_v10_0_5_18_3_3.json/g' \
  -e 's/release=10.0.5.18.3.2/release=10.0.5.18.3.3/g' \
  -e '/K_ramp_evaluated_inside_event_horizon=true/a\
single_trial_event_horizon_active=true\
state_dependent_log_rate_rejection_active=false\
PF_transport_operator_preserved=true\
target_aware_long_growth_corridor=true\
committed_target_propagated_before_mesh=true\
corridor_resolution_metric=h_tip_over_L_pz\
tip_h_over_da_role=audit_warning_only\
child_area_ratio_floor_relaxed=false' \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
