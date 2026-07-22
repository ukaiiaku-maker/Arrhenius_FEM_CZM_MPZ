#!/usr/bin/env bash
set -euo pipefail

# Corrected post-emission-refresh 1-D rerun with diagnostic work/energy output.
# MODE=smoke runs four spread-out registry rows; MODE=full runs all 96.

PYTHON_BIN="${PYTHON_BIN:-python}"
MODE="${MODE:-smoke}"
REGISTRY_96="${REGISTRY_96:-candidates/v9_12_1d_balanced_promotion_96_registry.csv}"
SMOKE_REGISTRY="${SMOKE_REGISTRY:-candidates/v9_12_corrected_energy_balanced96_smoke4_registry.csv}"
PROTOCOL="${PROTOCOL:-mpz_v9_12_protocol_example.csv}"
PHYSICS="${PHYSICS:-mpz_v9_12_emergent_gnd_common_physics.json}"
RESUME="${RESUME:-0}"

case "$MODE" in
  smoke)
    OUT="${OUT:-runs/v9_12_corrected_energy_balanced96_smoke4_c2_dt0p025_v1}"
    PROGRESS_EVERY="${PROGRESS_EVERY:-1}"
    ;;
  full)
    OUT="${OUT:-runs/v9_12_corrected_energy_balanced96_c2_dt0p025_v1}"
    PROGRESS_EVERY="${PROGRESS_EVERY:-8}"
    ;;
  *)
    echo "ERROR: MODE must be smoke or full; received: $MODE" >&2
    exit 1
    ;;
esac

if test ! -f "$REGISTRY_96"; then
  found="$({
    /usr/bin/find candidates -maxdepth 1 -type f \
      -name '*96*registry*.csv' -print 2>/dev/null || true
  } | /usr/bin/head -n 2)"
  found_count="$(printf '%s\n' "$found" | /usr/bin/awk 'NF {n += 1} END {print n + 0}')"
  if test "$found_count" = "1"; then
    REGISTRY_96="$found"
    echo "AUTO_SELECTED_REGISTRY $REGISTRY_96"
  else
    echo "ERROR: missing 96-case registry: $REGISTRY_96" >&2
    echo "Set REGISTRY_96 explicitly. Matching files found:" >&2
    printf '%s\n' "$found" >&2
    exit 1
  fi
fi

for required in "$REGISTRY_96" "$PROTOCOL" "$PHYSICS"; do
  test -e "$required" || {
    echo "ERROR: required input missing: $required" >&2
    exit 1
  }
done

"$PYTHON_BIN" -m py_compile \
  arrhenius_fracture/emergent_gnd_state_v912_energy.py \
  arrhenius_fracture/emergent_gnd_campaign_v912_energy.py \
  scripts/analyze_mpz_v9_12_energy_campaign.py \
  scripts/run_mpz_v9_12_emergent_gnd_screen_resilient.py

"$PYTHON_BIN" -m pytest -q \
  tests/test_emergent_gnd_dbtt_v912.py \
  tests/test_emergent_gnd_stiff_v912.py \
  tests/test_emergent_gnd_energy_v912.py

RUN_REGISTRY="$REGISTRY_96"
if test "$MODE" = "smoke"; then
  "$PYTHON_BIN" - "$REGISTRY_96" "$SMOKE_REGISTRY" <<'PY'
from pathlib import Path
import sys
import pandas as pd

source_path, out_path = sys.argv[1:]
frame = pd.read_csv(source_path, low_memory=False)
if len(frame) < 4:
    raise RuntimeError(f"need at least four candidates; found {len(frame)}")
indices = sorted(set([0, len(frame) // 3, 2 * len(frame) // 3, len(frame) - 1]))
if len(indices) != 4:
    raise RuntimeError(f"failed to choose four distinct smoke rows: {indices}")
selected = frame.iloc[indices].copy()
Path(out_path).parent.mkdir(parents=True, exist_ok=True)
selected.to_csv(out_path, index=False)
print(f"ENERGY_SMOKE_REGISTRY rows={len(selected)} out={out_path}")
print(selected["candidate_id"].to_string(index=False))
PY
  RUN_REGISTRY="$SMOKE_REGISTRY"
fi

if test "$RESUME" != "1"; then
  test ! -e "$OUT" || {
    echo "ERROR: output already exists: $OUT" >&2
    echo "Use RESUME=1 only for a compatible interrupted run, or choose a new OUT." >&2
    exit 1
  }
fi
mkdir -p "$OUT"

cmd=(
  "$PYTHON_BIN" -u
  scripts/run_mpz_v9_12_emergent_gnd_screen_resilient.py
  --stage 1d
  --candidate-registry "$RUN_REGISTRY"
  --protocol-csv "$PROTOCOL"
  --physics-json "$PHYSICS"
  --temperatures 300 400 500 600 700 800 900 1000 1100 1200
  --window-um 10 30
  --min-amplitude 8
  --target-localization 0.50
  --max-width-K 200
  --quiet-inner
  --progress-every "$PROGRESS_EVERY"
  --out "$OUT"
)
if test "$RESUME" = "1"; then
  cmd+=(--resume)
fi

echo "CORRECTED_ENERGY_CAMPAIGN_START mode=$MODE registry=$RUN_REGISTRY out=$OUT"
/usr/bin/caffeinate -dimsu \
/usr/bin/env \
PYTHONUNBUFFERED=1 \
MPZ_V912_COUPLED_OPERATOR_SUBSTEPS=2 \
MPZ_V912_MAX_FEEDBACK_SUBSTEP_S=0.025 \
"${cmd[@]}" \
2>&1 | /usr/bin/tee "$OUT/driver.log"

"$PYTHON_BIN" - "$OUT/resilient_progress.json" <<'PY'
from pathlib import Path
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
if payload.get("complete_count") != payload.get("total_candidates"):
    raise RuntimeError(f"campaign did not complete cleanly: {payload}")
if payload.get("unresolved_count") != 0:
    raise RuntimeError(f"campaign has unresolved candidates: {payload}")
print(
    "CORRECTED_ENERGY_COMPLETION_OK "
    f"candidates={payload['total_candidates']}"
)
PY

"$PYTHON_BIN" -u scripts/analyze_mpz_v9_12_energy_campaign.py \
  --root "$OUT" \
  --window-um 10 30 \
  --out-prefix "$OUT/energy_analysis"

echo "CORRECTED_ENERGY_CAMPAIGN_COMPLETE mode=$MODE out=$OUT"
