#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_2_four_class_focused_100um.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.2 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
export EMISSION_EVENT_HORIZON_FACTOR=${EMISSION_EVENT_HORIZON_FACTOR:-1.25}
export EMISSION_MAX_EVENTS_PER_HALF_STEP=${EMISSION_MAX_EVENTS_PER_HALF_STEP:-10000}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051832_single_trial_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission/mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_single_trial_stochastic_emission/g' \
  -e '/K_ramp_evaluated_inside_event_horizon=true/a\
single_trial_event_horizon_active=true\
state_dependent_log_rate_rejection_active=false\
log_rate_change_is_diagnostic_only=true\
PF_transport_operator_preserved=true\
events_batched=false\
event_horizon_factor='"$EMISSION_EVENT_HORIZON_FACTOR" \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
