#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_2_four_class_focused_100um.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.2 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
export EMISSION_INNER_MAX_ACTION_ERROR=${EMISSION_INNER_MAX_ACTION_ERROR:-0.01}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051832_action_error_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission/mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_action_error_stochastic_emission/g' \
  -e '/K_ramp_evaluated_inside_event_horizon=true/a\
action_error_control_active=true\
action_error_schema=v10.0.5.18.3.2_full_vs_two_half_integrated_hazard_action_error\
action_error_tolerance='"$EMISSION_INNER_MAX_ACTION_ERROR" \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
