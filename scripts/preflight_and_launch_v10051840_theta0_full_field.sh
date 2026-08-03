#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
EXPECTED_BRANCH=v10.0.5.18.4.0-theta0-pf-parity
PF_REFERENCE_ROOT=${PF_REFERENCE_ROOT:-/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1}
PARAMETER_SOURCE_ROOT=${PARAMETER_SOURCE_ROOT:-$ROOT/runtime_inputs/v10_0_5_18_four_class}
OPTION=v913_paper_peak01_0242980_persistent_sites
PF_CASE=${PF_CASE:-$PF_REFERENCE_ROOT/runs/v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/$OPTION/T1000K_th0_seed8666}
FAMILY_JSON=${FAMILY_JSON:-$PF_REFERENCE_ROOT/runs/v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json}
FAMILY_SHA_REQUIRED=a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a
RUN_ROOT=${RUN_ROOT:-$ROOT/runs/v10_0_5_18_4_0_peak1000_theta0_seed8666_pf_full_field_parity_v1}
PREFLIGHT_ROOT=${PREFLIGHT_ROOT:-$ROOT/runs/v10_0_5_18_4_0_theta0_full_field_preflight}
PREFLIGHT_OUT=$PREFLIGHT_ROOT/$OPTION/T1000K

cd "$ROOT"

BRANCH=$(git branch --show-current)
HEAD=$(git rev-parse HEAD)
if [[ "$BRANCH" != "$EXPECTED_BRANCH" ]]; then
  echo "ERROR: expected branch $EXPECTED_BRANCH; observed $BRANCH" >&2
  exit 1
fi
if [[ -n "$(git status --porcelain)" ]]; then
  echo "ERROR: working tree is not clean" >&2
  git status --short >&2
  exit 1
fi

for path in \
  "$PARAMETER_SOURCE_ROOT/v10_2_27_v913_four_class_paper_registry.csv" \
  "$PARAMETER_SOURCE_ROOT/four_class_transfer_manifest_v10_0_5_18.json" \
  "$FAMILY_JSON" \
  "$PF_CASE/steps_1000K.csv"
do
  [[ -f "$path" ]] || { echo "ERROR: missing required input: $path" >&2; exit 1; }
done

FAMILY_JSON=$(cd "$(dirname "$FAMILY_JSON")" && pwd)/$(basename "$FAMILY_JSON")
FAMILY_SHA=$(python - "$FAMILY_JSON" <<'PY'
import hashlib
import pathlib
import sys
print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)
if [[ "$FAMILY_SHA" != "$FAMILY_SHA_REQUIRED" ]]; then
  echo "ERROR: exact PF kernel SHA-256 mismatch" >&2
  echo "Expected: $FAMILY_SHA_REQUIRED" >&2
  echo "Observed: $FAMILY_SHA" >&2
  exit 1
fi

if [[ -e "$RUN_ROOT" ]]; then
  echo "ERROR: production run root already exists: $RUN_ROOT" >&2
  echo "Choose a new RUN_ROOT explicitly; the launcher will not overwrite a prior run." >&2
  exit 1
fi

python -m pip install -e '.[dev]'
python - <<'PY'
from importlib.metadata import version
observed = version('arrhenius-fem-czm')
assert observed == '10.0.5.18.4.0', observed
print('Package version verified:', observed)
PY

pytest -q \
  tests/test_theta0_pf_parity_v10051840.py \
  tests/test_theta0_pf_full_field_parity_v10051840.py \
  tests/test_theta0_pf_full_field_metadata_v10051840.py \
  tests/test_corridor_refinement_metadata_v10051840.py \
  tests/test_atomic_path_corridor_v10051839.py \
  tests/test_atomic_path_corridor_node_compaction_v10051839.py

rm -rf "$PREFLIGHT_ROOT"
mkdir -p "$PREFLIGHT_OUT"

export ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM=20
export ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY=0.035
export ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO=0.08
export ARRHENIUS_CORRIDOR_MAX_PATH_SUBSEGMENTS=16
export MIN_GLOBAL_FORWARD=0.05

export CLEAVAGE_HAZARD_MODE=exponential
export CLEAVAGE_HAZARD_SEED=8666
export CLEAVAGE_HAZARD_MIN_THRESHOLD=1e-12
export CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled
export CLEAVAGE_EVENT_MIN_FACTOR=0.5
export CLEAVAGE_EVENT_MAX_FACTOR=4.0
export EMISSION_HAZARD_SEED=8666
export EMISSION_HAZARD_MIN_THRESHOLD=1e-12
export EMISSION_MAX_ACTION_SUBSTEP=0.05
export EMISSION_MAX_EVENTS_PER_HALF_STEP=10000
export EMISSION_RAMP_MAX_LOG_RATE_CHANGE=0.25
export EMISSION_RAMP_ACTION_FLOOR=1e-6
export EMISSION_INNER_MAX_LOG_HAZARD_CHANGE=0.25
export EMISSION_INNER_ACTION_FLOOR=1e-6
export EMISSION_EVENT_HORIZON_FACTOR=1.25
export EMISSION_INNER_MAX_REFINEMENTS=24
export TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS=8

PARAMETER_SOURCE_ROOT="$PARAMETER_SOURCE_ROOT" \
FAMILY_JSON="$FAMILY_JSON" \
PREFLIGHT_OUT="$PREFLIGHT_OUT" \
python - <<'PY'
import json
import os
from pathlib import Path

from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production as production,
)
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity as control,
)

out = Path(os.environ['PREFLIGHT_OUT']).resolve()
production.main([
    '--parameter-source-root', os.environ['PARAMETER_SOURCE_ROOT'],
    '--parameter-option', control.REFERENCE_OPTION,
    '--signed-kernel-family', os.environ['FAMILY_JSON'],
    '--tip-refinement-radius-um', '330',
    '--selected-cluster-J-outer-um', '240',
    '--local-J-outer-um', '100',
    '--mode', '2d',
    '--bulk-plasticity-mode', 'full_field',
    '--temperatures', '1000',
    '--steps', '1',
    '--nx', '36', '--ny', '72',
    '--tip-h-fine', '1e-6', '--tip-ratio', '1.20',
    '--dU', '2e-7', '--dt', '8.4', '--n-stagger', '2',
    '--print-every', '1',
    '--adaptive-events',
    '--adaptive-event-target', '0.15',
    '--adaptive-min-frac', '1e-8',
    '--adaptive-grow', '4',
    '--da-phys', '5e-6',
    '--target-crack-extension-um', '20',
    '--crystal-aniso', '--crystal-compete',
    '--crystal-theta-deg', '0',
    '--min-global-forward', '0.05',
    '--crystal-C11', '523e9',
    '--crystal-C12', '203e9',
    '--crystal-C44', '160e9',
    '--cleave-gamma-aniso', '0.3',
    '--crystal-material', 'w',
    '--max-fronts', '1',
    '--crack-backend', 'adaptive_czm',
    '--czm-max-angle-error-deg', '35',
    '--j-decomposition', 'cluster',
    '--mpz-length-um', '50', '--mpz-n-bins', '80',
    '--save-snapshots', '1', '--snapshot-cols', '1',
    '--snapshot-by-crack-extension-um', '20',
    '--no-plots',
    '--out', str(out),
])

audit = json.loads((out / production.AUDIT_FILE).read_text())
assert audit['run_completed_without_exception'] is True, audit.get('runtime_error')
assert audit['bulk_parity']['FEM_bulk_plasticity_mode'] == 'full_field'
assert audit['bulk_parity']['detailed_balance_forward_minus_reverse'] is True

manifest = json.loads(
    (out / 'persistent_site_production_manifest_v10_0_5_18_3_9.json').read_text()
)
assert manifest['crystal_theta_deg'] == 0.0
assert manifest['bulk_state_evolves_in_fem'] is True
assert manifest['bulk_plasticity_semantic_mode'] == 'full_field'
lifecycle = manifest['refinement_metadata_lifecycle_audit']
assert lifecycle['captured_mesh_has_refinement_metadata'] is True
assert lifecycle['production_refinement_radius_m'] >= 330.0e-6

summary = json.loads((out / 'summary.json').read_text())[0]
assert summary['crystal_theta_deg'] == 0.0
assert summary['material_class'] == 'peak'
assert summary['bulk_state_evolves_in_fem'] is True
assert summary['taylor_renewal_time_s'] == 1.0e-9

metadata = json.loads((out / production.METADATA_AUDIT_FILE).read_text())
assert metadata['metadata_only'] is True
assert metadata['physics_or_mesh_modified'] is False
print('Real one-step theta=0 full-field production smoke passed.')
PY

mkdir -p "$RUN_ROOT"
cat > "$RUN_ROOT/local_preflight.json" <<EOF
{
  "schema": "local_preflight_v10_0_5_18_4_0",
  "git_branch": "$BRANCH",
  "git_head": "$HEAD",
  "parameter_option": "$OPTION",
  "temperature_K": 1000,
  "theta_deg": 0,
  "hazard_seed": 8666,
  "PF_case": "$PF_CASE",
  "PF_kernel": "$FAMILY_JSON",
  "PF_kernel_sha256": "$FAMILY_SHA",
  "focused_tests_passed": true,
  "real_one_step_full_field_smoke_passed": true,
  "physics_or_mesh_floors_relaxed": false
}
EOF

LOG=$RUN_ROOT/launcher.nohup.log
PIDFILE=$RUN_ROOT/launcher.pid
nohup env \
  ROOT="$ROOT" \
  PF_REFERENCE_ROOT="$PF_REFERENCE_ROOT" \
  PARAMETER_SOURCE_ROOT="$PARAMETER_SOURCE_ROOT" \
  PF_CASE="$PF_CASE" \
  FAMILY_JSON="$FAMILY_JSON" \
  CAMPAIGN_ROOT="$RUN_ROOT" \
  TARGET_EXT_UM=1000 \
  STEPS=2000000 \
  PRINT_EVERY=200 \
  SAVE_SNAPSHOTS=20 \
  SNAPSHOT_COLS=5 \
  SNAPSHOT_BY_EXT_UM=50 \
  GENERATE_SOLVER_PLOTS=1 \
  bash "$ROOT/run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh" \
  > "$LOG" 2>&1 &
PID=$!
echo "$PID" > "$PIDFILE"

sleep 5
if ! kill -0 "$PID" 2>/dev/null; then
  echo "ERROR: production launcher exited during startup" >&2
  tail -n 120 "$LOG" >&2 || true
  exit 1
fi

printf '\nProduction run launched successfully.\n'
echo "  Git head: $HEAD"
echo "  PID:      $PID"
echo "  Run root: $RUN_ROOT"
echo "  Log:      $LOG"
echo "  Case log: $RUN_ROOT/$OPTION/T1000K/console.log"
echo "  Analysis: $RUN_ROOT/parity_analysis"
