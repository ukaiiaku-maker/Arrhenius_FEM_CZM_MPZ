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
export ARRHENIUS_EVENT_RESOLVED_PATCH_LEVELS=${ARRHENIUS_EVENT_RESOLVED_PATCH_LEVELS:-12}
export ARRHENIUS_EVENT_RESOLVED_MAX_RAY_SEGMENTS=${ARRHENIUS_EVENT_RESOLVED_MAX_RAY_SEGMENTS:-512}
export CLEAVAGE_EVENT_MIN_FACTOR=${CLEAVAGE_EVENT_MIN_FACTOR:-0.5}
export CLEAVAGE_EVENT_MAX_FACTOR=${CLEAVAGE_EVENT_MAX_FACTOR:-4.0}

TMP=$(mktemp "${TMPDIR:-/tmp}/v10051838_centroid_cavity_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope/mode_i_first_passage_v10_0_5_18_3_8_centroid_cavity/g' \
  -e 's/persistent_site_production_manifest_v10_0_5_18_3_4.json/persistent_site_production_manifest_v10_0_5_18_3_8.json/g' \
  -e 's/release=10.0.5.18.3.4/release=10.0.5.18.3.8/g' \
  -e 's/Focused v10.0.5.18.3.4 path-envelope validation complete/Focused v10.0.5.18.3.8 centroid-cavity validation complete/g' \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
