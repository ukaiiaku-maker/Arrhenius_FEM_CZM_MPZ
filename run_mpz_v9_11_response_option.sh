#!/usr/bin/env bash
set -euo pipefail

# Named response options preserve scientifically distinct mechanisms instead of
# treating every DBTT-like curve as evidence for the same physics.  The default
# is the production DBTT candidate; peak, broad-shielding, intrinsic-control,
# moderate-shielding, and weak-temperature options are selected with OPTION.
OPTION=${OPTION:-dbtt_primary}
CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm/bin/python}
CANONICAL_ROOT=${CANONICAL_ROOT:-mpz_v9_11_parameters}
OPTION_ROOT_BASE=${OPTION_ROOT_BASE:-mpz_v9_11_response_option_roots}
PARAMETER_ROOT=${PARAMETER_ROOT:-$OPTION_ROOT_BASE/$OPTION}
BULK_PLASTICITY_MODE=${BULK_PLASTICITY_MODE:-tip_only}
TARGET_EXT_UM=${TARGET_EXT_UM:-500}
STEPS=${STEPS:-25000}
HEARTBEAT_SECONDS=${HEARTBEAT_SECONDS:-60}
TEMPS=${TEMPS:-}
OUTROOT=${OUTROOT:-runs/mpz_v9_11_1_${OPTION}_${BULK_PLASTICITY_MODE}_Rcurve_v1}

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: Python executable not found: $PYTHON_BIN" >&2
  exit 2
fi

OPTION_FIELDS=$(
  "$PYTHON_BIN" - "$OPTION" <<'PY'
import json
import sys
from pathlib import Path

option = sys.argv[1]
registry = json.loads(Path("mpz_v9_11_response_options.json").read_text())
if option not in registry["options"]:
    raise SystemExit(f"unknown response option: {option}")
row = registry["options"][option]
fields = [
    row["material_class"],
    row["candidate_id"],
    row["role"],
    " ".join(str(x) for x in row.get("recommended_followup_temperatures_K", [])),
]
print("\t".join(fields))
PY
)
IFS=$'\t' read -r CLASS CANDIDATE_ID OPTION_ROLE RECOMMENDED_TEMPS <<< "$OPTION_FIELDS"

if [[ -z "$TEMPS" ]]; then
  if [[ -n "$RECOMMENDED_TEMPS" ]]; then
    TEMPS=$RECOMMENDED_TEMPS
  else
    TEMPS="300 700 900 1200"
  fi
fi

if [[ ! -f "$PARAMETER_ROOT/response_option_selection.json" ]]; then
  "$PYTHON_BIN" prepare_mpz_v9_11_response_option.py \
    --option "$OPTION" \
    --canonical-root "$CANONICAL_ROOT" \
    --outroot "$PARAMETER_ROOT"
fi

mkdir -p "$OUTROOT"
REPORT_LOG=$OUTROOT/reported_campaign.log

printf '[%s] CAMPAIGN_START option=%s candidate=%s role=%s class=%s temps="%s" target=%sum steps=%s outroot=%s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$OPTION" "$CANDIDATE_ID" "$OPTION_ROLE" \
  "$CLASS" "$TEMPS" "$TARGET_EXT_UM" "$STEPS" "$OUTROOT" | tee "$REPORT_LOG"

set +e
PYTHONUNBUFFERED=1 \
CONDA_ENV="$CONDA_ENV" \
PYTHON_BIN="$PYTHON_BIN" \
PARAMETER_ROOT="$PARAMETER_ROOT" \
CLASS="$CLASS" \
BULK_PLASTICITY_MODE="$BULK_PLASTICITY_MODE" \
TEMPS="$TEMPS" \
OUTROOT="$OUTROOT" \
TARGET_EXT_UM="$TARGET_EXT_UM" \
STEPS="$STEPS" \
MPZ_LENGTH_UM="${MPZ_LENGTH_UM:-50}" \
MPZ_N_BINS="${MPZ_N_BINS:-80}" \
SKIP_EXISTING="${SKIP_EXISTING:-1}" \
bash run_mpz_v9_11_mode_i_rcurve_3T.sh \
  >> "$REPORT_LOG" 2>&1 &
runner_pid=$!
set -e

while kill -0 "$runner_pid" 2>/dev/null; do
  active_log=$(find "$OUTROOT" -type f -name run.log -print0 2>/dev/null \
    | xargs -0 ls -1t 2>/dev/null | head -n 1 || true)
  if [[ -n "$active_log" ]]; then
    progress=$(tail -n 250 "$active_log" 2>/dev/null \
      | grep '\[T=' | tail -n 1 || true)
    if [[ -n "$progress" ]]; then
      printf '[%s] HEARTBEAT option=%s pid=%s %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$OPTION" "$runner_pid" "$progress"
    else
      printf '[%s] HEARTBEAT option=%s pid=%s active_log=%s waiting_for_progress\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$OPTION" "$runner_pid" "$active_log"
    fi
  else
    printf '[%s] HEARTBEAT option=%s pid=%s initializing\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$OPTION" "$runner_pid"
  fi
  sleep "$HEARTBEAT_SECONDS"
done

set +e
wait "$runner_pid"
rc=$?
set -e

summary=$OUTROOT/rcurve_temperature_summary.csv
if [[ "$rc" -eq 0 ]]; then
  printf '[%s] CAMPAIGN_COMPLETE option=%s returncode=0 summary=%s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S')" "$OPTION" "$summary" | tee -a "$REPORT_LOG"
else
  printf '[%s] CAMPAIGN_INCOMPLETE_OR_FAILED option=%s returncode=%s summary=%s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S')" "$OPTION" "$rc" "$summary" | tee -a "$REPORT_LOG" >&2
fi

if [[ -f "$summary" ]]; then
  "$PYTHON_BIN" - "$summary" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
print("FINAL_CASE_STATUS")
with path.open(newline="") as handle:
    for row in csv.DictReader(handle):
        print(
            f"T={row.get('T_K')}K status={row.get('status')} "
            f"extension={row.get('final_extension_um')}um "
            f"events={row.get('n_growth_events')} "
            f"Kinit={row.get('K_init_MPa_sqrt_m')}"
        )
PY
fi

exit "$rc"
