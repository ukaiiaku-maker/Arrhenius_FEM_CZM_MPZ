#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN=${PYTHON_BIN:-python}
MODE=${MODE:-smoke}
REGISTRY=${REGISTRY:-candidates/v9_13_persistent_sites_top5_registry.csv}
PROTOCOL=${PROTOCOL:-mpz_v9_12_protocol_example.csv}
BASE_PHYSICS=${BASE_PHYSICS:-mpz_v9_13_persistent_sites_common_physics.json}
SIGNED_KERNEL_FAMILY_JSON=${SIGNED_KERNEL_FAMILY_JSON:-}
ALLOW_UNIT_LINE_CONVERSION=${ALLOW_UNIT_LINE_CONVERSION:-0}
PHYSICAL_MIN_FRONT_WIDTH_NM=${PHYSICAL_MIN_FRONT_WIDTH_NM:-}
COUPLED_MOVING_TIP_AUDIT=${COUPLED_MOVING_TIP_AUDIT:-0}

case "$COUPLED_MOVING_TIP_AUDIT" in
  0|1) ;;
  *)
    echo "ERROR: COUPLED_MOVING_TIP_AUDIT must be 0 or 1" >&2
    exit 2
    ;;
esac

case "$MODE" in
  smoke)
    TEMPERATURES=${TEMPERATURES:-"700 800 900 1000 1100 1200"}
    ROWS=1
    SUBSTEPS=${SUBSTEPS:-2}
    OUT=${OUT:-runs/v9_13_persistent_sites_top1_smoke_c2_dt0p025_v1}
    ;;
  full)
    TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 900 1000 1100 1200"}
    ROWS=5
    SUBSTEPS=${SUBSTEPS:-2}
    OUT=${OUT:-runs/v9_13_persistent_sites_top5_c2_dt0p025_v1}
    ;;
  convergence)
    TEMPERATURES=${TEMPERATURES:-"700 800 900 1000 1100 1200"}
    ROWS=2
    SUBSTEPS=${SUBSTEPS:-4}
    OUT=${OUT:-runs/v9_13_persistent_sites_top2_c4_dt0p025_v1}
    ;;
  *)
    echo "ERROR: MODE must be smoke, full, or convergence" >&2
    exit 2
    ;;
esac

for required in "$REGISTRY" "$PROTOCOL" "$BASE_PHYSICS"; do
  test -f "$required" || { echo "ERROR: missing $required" >&2; exit 2; }
done

test ! -e "$OUT" || {
  echo "ERROR: output already exists: $OUT" >&2
  exit 2
}
mkdir -p "$OUT/runtime_inputs"

RUN_PHYSICS="$BASE_PHYSICS"
if test -n "$SIGNED_KERNEL_FAMILY_JSON"; then
  test -f "$SIGNED_KERNEL_FAMILY_JSON" || {
    echo "ERROR: signed kernel family not found: $SIGNED_KERNEL_FAMILY_JSON" >&2
    exit 2
  }
  RUN_PHYSICS="$OUT/runtime_inputs/v9_13_physics_with_2d_line_conversion.json"
  "$PYTHON_BIN" scripts/extract_v10221_line_conversion_for_v913.py \
    --signed-kernel-family "$SIGNED_KERNEL_FAMILY_JSON" \
    --base-physics "$BASE_PHYSICS" \
    --out "$RUN_PHYSICS"
elif test "$ALLOW_UNIT_LINE_CONVERSION" != 1; then
  cat >&2 <<'EOF'
ERROR: exact 1-D/2-D matching requires the signed 2-D kernel family so the
activation_to_line_content normalization can be imported.

Set, for example:
  SIGNED_KERNEL_FAMILY_JSON=/path/to/v10_2_14_active_only_campaign_family.json

For a source-law-only diagnostic using the legacy 1-D unit line conversion,
set ALLOW_UNIT_LINE_CONVERSION=1 explicitly.
EOF
  exit 2
else
  echo "WARNING: using legacy 1-D unit activation-to-line conversion" >&2
fi

if test -n "$PHYSICAL_MIN_FRONT_WIDTH_NM" \
  || test "$COUPLED_MOVING_TIP_AUDIT" = 1; then
  CONFIGURED_PHYSICS="$OUT/runtime_inputs/v9_13_configured_physics.json"
  "$PYTHON_BIN" - \
    "$RUN_PHYSICS" \
    "$CONFIGURED_PHYSICS" \
    "$PHYSICAL_MIN_FRONT_WIDTH_NM" \
    "$COUPLED_MOVING_TIP_AUDIT" <<'PY'
import json
from pathlib import Path
import sys

source, target, width_nm, coupled = sys.argv[1:]
payload = json.loads(Path(source).read_text())
common = payload.setdefault("common_physics", {})
if width_nm:
    width_m = float(width_nm) * 1.0e-9
    if not width_m > 0.0:
        raise ValueError("PHYSICAL_MIN_FRONT_WIDTH_NM must be positive")
    common["minimum_front_width_m"] = width_m
    payload["front_width_provenance"] = (
        f"explicit_physical_sensitivity_{float(width_nm):g}_nm"
    )
if coupled not in {"0", "1"}:
    raise ValueError("COUPLED_MOVING_TIP_AUDIT must be 0 or 1")
common["coupled_moving_tip_enabled"] = coupled == "1"
Path(target).write_text(json.dumps(payload, indent=2) + "\n")
print(
    "V913_CONFIGURED_PHYSICS "
    f"minimum_front_width_m={common.get('minimum_front_width_m')} "
    f"coupled_moving_tip_enabled={common['coupled_moving_tip_enabled']} "
    f"out={target}"
)
PY
  RUN_PHYSICS="$CONFIGURED_PHYSICS"
fi

RUN_REGISTRY="$OUT/runtime_inputs/selected_registry.csv"
"$PYTHON_BIN" - "$REGISTRY" "$RUN_REGISTRY" "$ROWS" <<'PY'
import csv
from pathlib import Path
import sys
source, out, count = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(source, newline="") as stream:
    reader = csv.DictReader(stream)
    rows = list(reader)
    fields = list(reader.fieldnames or [])
if len(rows) < count:
    raise RuntimeError(f"registry has {len(rows)} rows but {count} are required")
Path(out).parent.mkdir(parents=True, exist_ok=True)
with open(out, "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows[:count])
print(f"V913_SELECTED_REGISTRY rows={count} out={out}")
PY

"$PYTHON_BIN" -m py_compile \
  arrhenius_fracture/emergent_gnd_types_v913.py \
  arrhenius_fracture/emergent_gnd_state_v913.py \
  arrhenius_fracture/emergent_gnd_campaign_v913.py \
  arrhenius_fracture/emergent_gnd_dbtt_v913.py \
  scripts/extract_v10221_line_conversion_for_v913.py \
  scripts/run_mpz_v9_13_persistent_top5.py

"$PYTHON_BIN" -m pytest -q \
  tests/test_emergent_gnd_persistent_v913.py \
  tests/test_emergent_gnd_stiff_v912.py \
  tests/test_emergent_gnd_energy_v912.py

MPZ_V912_COUPLED_OPERATOR_SUBSTEPS="$SUBSTEPS" \
MPZ_V912_MAX_FEEDBACK_SUBSTEP_S=0.025 \
/usr/bin/caffeinate -dimsu \
"$PYTHON_BIN" -u scripts/run_mpz_v9_13_persistent_top5.py \
  --candidate-registry "$RUN_REGISTRY" \
  --protocol-csv "$PROTOCOL" \
  --physics-json "$RUN_PHYSICS" \
  --temperatures $TEMPERATURES \
  --window-um 10 30 \
  --min-amplitude 8 \
  --target-localization 0.50 \
  --max-width-K 200 \
  --out "$OUT" \
  2>&1 | /usr/bin/tee "$OUT/driver.log"

echo "V913_PERSISTENT_TOP5_RUN_COMPLETE mode=$MODE out=$OUT"
