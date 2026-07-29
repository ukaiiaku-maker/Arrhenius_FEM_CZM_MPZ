#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_2_four_class_focused_100um.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.2 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
export MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}
export ARRHENIUS_CORRIDOR_GUARD_UM=${ARRHENIUS_CORRIDOR_GUARD_UM:-10}
export ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ=${ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ:-0.25}
export ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ=${ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ:-0.245}
export ARRHENIUS_ENVELOPE_BUFFER_LPZ=${ARRHENIUS_ENVELOPE_BUFFER_LPZ:-1.0}
export ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES=${ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES:-181}
export ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES=${ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES:-4}
export ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR=${ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR:-0.9}
export ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ=${ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ:-0.5}
export ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS=${ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS:-16}
export EMISSION_EVENT_HORIZON_FACTOR=${EMISSION_EVENT_HORIZON_FACTOR:-1.25}
export EMISSION_MAX_EVENTS_PER_HALF_STEP=${EMISSION_MAX_EVENTS_PER_HALF_STEP:-10000}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051834_path_envelope_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission/mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope/g' \
  -e 's/persistent_site_production_manifest_v10_0_5_18_3_2.json/persistent_site_production_manifest_v10_0_5_18_3_4.json/g' \
  -e 's/release=10.0.5.18.3.2/release=10.0.5.18.3.4/g' \
  -e '/theta_deg=\$THETA/a\
minimum_global_forward_cosine=$MIN_GLOBAL_FORWARD' \
  -e '/--crystal-theta-deg "\$THETA"/a\
      --min-global-forward "$MIN_GLOBAL_FORWARD" \\' \
  -e '/K_ramp_evaluated_inside_event_horizon=true/a\
single_trial_event_horizon_active=true\
state_dependent_log_rate_rejection_active=false\
PF_transport_operator_preserved=true\
path_aware_growth_envelope=true\
directional_support_basis=complete_global_forward_gate_conservative\
full_reachable_envelope_process_zone_resolved=true\
envelope_resolution_metric=nearest_fixed_element_patch_mean_edge_over_L_pz\
tip_h_over_da_role=audit_warning_only\
initial_physical_refinement_radius_preserved=true\
swept_physical_refinement_active=false\
child_area_ratio_floor_relaxed=false' \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
