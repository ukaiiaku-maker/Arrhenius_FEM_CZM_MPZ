#!/usr/bin/env bash
set -euo pipefail

# Always run the source tree that contains this script.  This prevents an
# older installed/adjacent arrhenius_fracture package from silently winning.
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1

# Remove stale bytecode from earlier overlays.
find arrhenius_fracture tests -type d -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null || true

python -B - <<'PY'
from pathlib import Path
import arrhenius_fracture.sn_pd2d_stateful as driver
import arrhenius_fracture.stateful_peridynamics as pdmod

root = Path.cwd().resolve()
driver_path = Path(driver.__file__).resolve()
pd_path = Path(pdmod.__file__).resolve()
expected = "SN_2D_intact_FEM_stateful_local_peridynamics_v3_1_root_localized"

print(f"STATEFUL_PD preflight project_root: {root}")
print(f"STATEFUL_PD preflight driver:       {driver_path}")
print(f"STATEFUL_PD preflight PD module:    {pd_path}")
print(f"STATEFUL_PD preflight model:        {getattr(driver, 'MODEL_ID', None)}")

if root not in driver_path.parents or root not in pd_path.parents:
    raise SystemExit("ERROR: Python imported stateful-PD code outside PROJECT_ROOT")
if getattr(driver, "MODEL_ID", None) != expected:
    raise SystemExit("ERROR: active driver is not the v3.1 root-localized model")

dests = {a.dest for a in driver.build_parser()._actions}
required = {
    "pd_initiation_radius_m",
    "pd_initiation_taper_m",
    "pd_initiation_back_extent_m",
    "pd_amplification_damage_scale",
}
missing = sorted(required - dests)
if missing:
    raise SystemExit(f"ERROR: active parser is missing v3 options: {missing}")
print("STATEFUL_PD v3.1 preflight OK")
PY

OUT="${OUT:-runs/sn_stateful_pd_pilot}"
python -B -m arrhenius_fracture.sn_pd2d_stateful \
  --out "$OUT" \
  --cases ${CASES:-no_shield shielded} \
  --T "${T:-300}" --R "${R:-0.1}" --frequency-Hz "${FREQ:-1000}" \
  --sigma-a-MPa ${STRESSES:-500 600 700} \
  --cycles-max "${CYCLES_MAX:-1e9}" \
  --block-cycles "${BLOCK_CYCLES:-1e7}" \
  --max-blocks "${MAX_BLOCKS:-3000}" \
  --target-dep-eq-block "${TARGET_DEP_BLOCK:-2e-4}" \
  --target-rho-rel-block "${TARGET_RHO_REL:-0.05}" \
  --max-transition-probability "${MAX_TRANSITION_P:-0.08}" \
  --nx "${NX:-36}" --ny "${NY:-72}" \
  --root-h-fine "${ROOT_H_FINE:-30e-6}" \
  --pd-horizon-m "${PD_HORIZON:-90e-6}" \
  --pd-patch-radius-m "${PD_PATCH_RADIUS:-0.45e-3}" \
  --pd-boundary-shell-m "${PD_BOUNDARY_SHELL:-100e-6}" \
  --pd-initiation-radius-m "${PD_INITIATION_RADIUS:-240e-6}" \
  --pd-initiation-taper-m "${PD_INITIATION_TAPER:-60e-6}" \
  --pd-initiation-back-extent-m "${PD_INITIATION_BACK_EXTENT:-60e-6}" \
  --pd-amplification-damage-scale "${PD_AMPLIFICATION_DAMAGE_SCALE:-0.05}" \
  --site-density-m2 "${SITE_DENSITY_M2:-5e10}" \
  --hit-count "${HIT_COUNT:-3}" \
  --hit-memory-s "${HIT_MEMORY_S:-1e-6}" \
  --birth-scale "${BIRTH_SCALE:-1.0}" \
  --seed "${MESH_SEED:-1}" \
  --pd-seed "${PD_SEED:-${MESH_SEED:-1}}" \
  --softening-damage "${SOFTENING_DAMAGE:-1e-3}" \
  --plastic-n-phase "${PLASTIC_PHASES:-12}" \
  --hazard-n-phase "${HAZARD_PHASES:-16}" \
  --snapshot-every "${SNAPSHOT_EVERY:-25}" \
  --print-every "${PRINT_EVERY:-1}"
