#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_4_four_class_focused_long_growth.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.4 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
export ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY=${ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY:-0.035}
export ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO=${ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO:-0.08}
export CLEAVAGE_EVENT_MIN_FACTOR=${CLEAVAGE_EVENT_MIN_FACTOR:-0.5}
export CLEAVAGE_EVENT_MAX_FACTOR=${CLEAVAGE_EVENT_MAX_FACTOR:-4.0}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051836_committed_tip_audit_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope/mode_i_first_passage_v10_0_5_18_3_6_committed_tip_resolution_audit/g' \
  -e 's/persistent_site_production_manifest_v10_0_5_18_3_4.json/persistent_site_production_manifest_v10_0_5_18_3_6.json/g' \
  -e 's/release=10.0.5.18.3.4/release=10.0.5.18.3.6/g' \
  -e 's/Focused v10.0.5.18.3.4 path-envelope validation complete/Focused v10.0.5.18.3.6 committed-tip audit complete/g' \
  -e '/child_area_ratio_floor_relaxed=false/a\
committed_tip_resolution_audit_active=true\
committed_tip_and_candidate_endpoint_metrics_separated=true\
joint_triangle_quality_and_child_area_margins_recorded=true\
event_length_clip_state_recorded=true\
tip_h_fine_contract_recorded=true\
telemetry_only=true\
geometry_changed=false\
hazard_or_event_length_changed=false' \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
