#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_4_four_class_focused_long_growth.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.4 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
export ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS=${ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS:-8}
export ARRHENIUS_QUALITY_AWARE_PARTITIONS=${ARRHENIUS_QUALITY_AWARE_PARTITIONS:-2,4,8}
export ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY=${ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY:-0.035}
export ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO=${ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO:-0.08}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051835_quality_aware_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope/mode_i_first_passage_v10_0_5_18_3_5_four_class_quality_aware_czm_retry/g' \
  -e 's/persistent_site_production_manifest_v10_0_5_18_3_4.json/persistent_site_production_manifest_v10_0_5_18_3_5.json/g' \
  -e 's/release=10.0.5.18.3.4/release=10.0.5.18.3.5/g' \
  -e 's/Focused v10.0.5.18.3.4 path-envelope validation complete/Focused v10.0.5.18.3.5 quality-aware CZM validation complete/g' \
  -e '/child_area_ratio_floor_relaxed=false/a\
quality_gate_inside_retry_transaction=true\
shape_regular_local_tip_patch_bisection=true\
immediate_parent_child_area_ratio_enforced=true\
cumulative_child_area_ratio_role=audit_only\
quality_aware_patch_levels=$ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS\
quality_aware_partition_counts=$ARRHENIUS_QUALITY_AWARE_PARTITIONS' \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
