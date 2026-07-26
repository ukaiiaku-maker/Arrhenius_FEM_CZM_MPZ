#!/usr/bin/env bash
set -uo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_18_2_diagnostic_smoke}
LAUNCH_LOG=${LAUNCH_LOG:-$CAMPAIGN_ROOT/launcher_console.log}
STATUS_JSON=${STATUS_JSON:-$CAMPAIGN_ROOT/launcher_exit_status.json}
RUNNER="$ROOT/run_v10_0_5_18_2_four_class_focused_100um.sh"

mkdir -p "$CAMPAIGN_ROOT"
START_UTC=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

printf '[DIAGNOSTIC_START] utc=%s pid=%s runner=%s\n' \
  "$START_UTC" "$$" "$RUNNER" | tee "$LAUNCH_LOG"

set +e
PYTHONUNBUFFERED=1 bash "$RUNNER" 2>&1 | tee -a "$LAUNCH_LOG"
RUNNER_RC=${PIPESTATUS[0]}
set -e

END_UTC=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
SIGNAL=0
if (( RUNNER_RC >= 128 )); then
  SIGNAL=$((RUNNER_RC - 128))
fi

RUNNER_RC="$RUNNER_RC" \
SIGNAL="$SIGNAL" \
START_UTC="$START_UTC" \
END_UTC="$END_UTC" \
LAUNCH_LOG="$LAUNCH_LOG" \
STATUS_JSON="$STATUS_JSON" \
python - <<'PY'
import json
import os
from pathlib import Path

payload = {
    "schema": "v10.0.5.18.2_diagnostic_launcher_exit_status",
    "runner_exit_code": int(os.environ["RUNNER_RC"]),
    "terminating_signal": int(os.environ["SIGNAL"]),
    "started_utc": os.environ["START_UTC"],
    "ended_utc": os.environ["END_UTC"],
    "launcher_log": os.environ["LAUNCH_LOG"],
    "clean_exit": int(os.environ["RUNNER_RC"]) == 0,
}
Path(os.environ["STATUS_JSON"]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(payload, indent=2, sort_keys=True))
PY

if (( RUNNER_RC != 0 )); then
  printf '[DIAGNOSTIC_FAILURE] runner_exit_code=%s signal=%s\n' \
    "$RUNNER_RC" "$SIGNAL" | tee -a "$LAUNCH_LOG" >&2
fi

exit "$RUNNER_RC"
