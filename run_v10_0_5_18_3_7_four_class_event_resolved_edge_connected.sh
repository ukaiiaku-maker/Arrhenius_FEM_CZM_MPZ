#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BASE_RUNNER="$ROOT/run_v10_0_5_18_3_7_four_class_event_resolved_local_cavity.sh"

if [[ ! -f "$BASE_RUNNER" ]]; then
  echo "ERROR: missing base v10.0.5.18.3.7 runner: $BASE_RUNNER" >&2
  exit 1
fi

export ROOT
TMP=$(mktemp "${TMPDIR:-/tmp}/v10051837_edge_connected_runner.XXXXXX")
trap 'rm -f "$TMP"' EXIT

sed \
  -e 's/mode_i_first_passage_v10_0_5_18_3_7_event_resolved_local_cavity/mode_i_first_passage_v10_0_5_18_3_7_event_resolved_edge_connected/g' \
  -e 's/Focused v10.0.5.18.3.7 event-resolved local-cavity validation complete/Focused v10.0.5.18.3.7 edge-connected exact-ray validation complete/g' \
  "$BASE_RUNNER" > "$TMP"

chmod +x "$TMP"
bash "$TMP"
