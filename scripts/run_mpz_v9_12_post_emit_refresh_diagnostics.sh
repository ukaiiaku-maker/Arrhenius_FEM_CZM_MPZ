#!/usr/bin/env bash
set -euo pipefail

# Six-candidate convergence matrix for the v9.12 post-emission rate refresh.
# Run from the repository root in arrhenius-fem-czm-v912-gnd.

PYTHON_BIN="${PYTHON_BIN:-python}"
SOURCE_REGISTRY="${SOURCE_REGISTRY:-candidates/v9_12_targeted_finalists_11_registry.csv}"
DIAG_REGISTRY="${DIAG_REGISTRY:-candidates/v9_12_post_emit_refresh_diag6_registry.csv}"
PROTOCOL="${PROTOCOL:-mpz_v9_12_protocol_example.csv}"
PHYSICS="${PHYSICS:-mpz_v9_12_emergent_gnd_common_physics.json}"
OUT_PREFIX="${OUT_PREFIX:-runs/v9_12_post_emit_refresh_diag6}"
EXPECTED_INTEGRATOR="coupled_mobile_retained_backward_euler_v2_post_emit_refresh"

TEMPERATURES=(300 400 500 600 700 800 900 1000 1100 1200)

for required in "$SOURCE_REGISTRY" "$PROTOCOL" "$PHYSICS"; do
  test -e "$required" || {
    echo "ERROR: required input missing: $required" >&2
    exit 1
  }
done

"$PYTHON_BIN" -m py_compile \
  arrhenius_fracture/emergent_gnd_state_v912_stiff.py \
  scripts/run_mpz_v9_12_emergent_gnd_screen_resilient.py

"$PYTHON_BIN" -m pytest -q \
  tests/test_emergent_gnd_dbtt_v912.py \
  tests/test_emergent_gnd_stiff_v912.py

"$PYTHON_BIN" - "$SOURCE_REGISTRY" "$DIAG_REGISTRY" <<'PY'
from pathlib import Path
import sys
import pandas as pd

source_path, out_path = sys.argv[1:]
ids = [
    "v912_targeted_local_peak_013476_0314",
    "v912_targeted_local_peak_013476_0368",
    "v912_targeted_local_peak_013476_0396",
    "v912_targeted_local_peak_005518_0118",
    "v912_targeted_local_peak_013476_0274",
    "v912_targeted_local_plateau_010759_0403",
]

frame = pd.read_csv(source_path, low_memory=False)
frame["candidate_id"] = frame["candidate_id"].astype(str)
order = {candidate_id: index for index, candidate_id in enumerate(ids)}
selected = frame[frame["candidate_id"].isin(order)].copy()
missing = sorted(set(ids) - set(selected["candidate_id"]))
if missing:
    raise RuntimeError(f"missing diagnostic candidates: {missing}")
selected["_order"] = selected["candidate_id"].map(order)
selected = selected.sort_values("_order").drop(columns="_order")
if len(selected) != len(ids) or not selected["candidate_id"].is_unique:
    raise RuntimeError("diagnostic registry is incomplete or nonunique")
Path(out_path).parent.mkdir(parents=True, exist_ok=True)
selected.to_csv(out_path, index=False)
print(f"DIAGNOSTIC_REGISTRY rows={len(selected)} out={out_path}")
PY

run_setting () {
  local tag="$1"
  local coupled_substeps="$2"
  local feedback_dt="$3"
  local out="${OUT_PREFIX}_${tag}_v1"

  test ! -e "$out" || {
    echo "ERROR: output already exists: $out" >&2
    return 1
  }
  mkdir -p "$out"

  echo \
    "POST_EMIT_DIAGNOSTIC_START tag=$tag coupled_substeps=$coupled_substeps "\
    "feedback_dt=$feedback_dt out=$out"

  /usr/bin/caffeinate -dimsu \
  /usr/bin/env \
  PYTHONUNBUFFERED=1 \
  MPZ_V912_COUPLED_OPERATOR_SUBSTEPS="$coupled_substeps" \
  MPZ_V912_MAX_FEEDBACK_SUBSTEP_S="$feedback_dt" \
  "$PYTHON_BIN" -u \
    scripts/run_mpz_v9_12_emergent_gnd_screen_resilient.py \
    --stage 1d \
    --candidate-registry "$DIAG_REGISTRY" \
    --protocol-csv "$PROTOCOL" \
    --physics-json "$PHYSICS" \
    --temperatures "${TEMPERATURES[@]}" \
    --window-um 10 30 \
    --min-amplitude 50 \
    --target-localization 0.50 \
    --max-width-K 200 \
    --quiet-inner \
    --progress-every 1 \
    --out "$out" \
    2>&1 | /usr/bin/tee "$out/campaign.log"

  "$PYTHON_BIN" - "$out" "$EXPECTED_INTEGRATOR" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
expected = sys.argv[2]
summaries = sorted(root.glob("*/candidate_summary.json"))
if len(summaries) != 6:
    raise RuntimeError(f"expected 6 summaries in {root}; found {len(summaries)}")
for path in summaries:
    payload = json.loads(path.read_text())
    metadata = payload.get("integration_metadata", {})
    observed = metadata.get("spatial_integrator")
    if observed != expected:
        raise RuntimeError(
            f"wrong integrator metadata in {path}: {observed!r} != {expected!r}"
        )
    if metadata.get("constitutive_feedback_update") != (
        "refresh_after_midpoint_emission"
    ):
        raise RuntimeError(f"missing post-emission refresh metadata in {path}")
print(f"INTEGRATOR_METADATA_OK summaries={len(summaries)} root={root}")
PY
}

run_setting c2_dt0p1_postemit_v2   2 0.1
run_setting c2_dt0p05_postemit_v2  2 0.05
run_setting c2_dt0p025_postemit_v2 2 0.025
run_setting c4_dt0p05_postemit_v2  4 0.05

echo "POST_EMIT_DIAGNOSTICS_COMPLETE"
