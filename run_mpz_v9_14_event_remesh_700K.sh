#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT=${OUTROOT:-runs/mpz_v9_14_event_remesh_700K_v1}
SEEDS=${SEEDS:-"1"}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
T_K=${T_K:-700}
TARGET_EXT_UM=${TARGET_EXT_UM:-50}
STEPS=${STEPS:-12000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
N_STAGGER=${N_STAGGER:-2}
PRINT_EVERY=${PRINT_EVERY:-25}
# v9.14 interprets this as an absolute integrated-hazard increment dB,
# not a fraction of the remaining renewal threshold.
ADAPTIVE_EVENT_TARGET=${ADAPTIVE_EVENT_TARGET:-0.01}
DA_PHYS_UM=${DA_PHYS_UM:-5}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-45}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-5}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-10}
EVENT_STATISTICS=${EVENT_STATISTICS:-deterministic}
STOCHASTIC_EMISSION=${STOCHASTIC_EMISSION:-0}
EVENT_REMESH_TARGET_H_M=${EVENT_REMESH_TARGET_H_M:-1e-6}
EVENT_REMESH_PATCH_RADIUS_UM=${EVENT_REMESH_PATCH_RADIUS_UM:-25}
EVENT_REMESH_MAX_EDGE_SPLITS=${EVENT_REMESH_MAX_EDGE_SPLITS:-256}
EVENT_REMESH_TARGET_EDGE_FACTOR=${EVENT_REMESH_TARGET_EDGE_FACTOR:-1.25}
EVENT_REMESH_MIN_QUALITY=${EVENT_REMESH_MIN_QUALITY:-0.02}
SKIP_EXISTING=${SKIP_EXISTING:-1}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: Python executable not found: $PYTHON_BIN" >&2
  exit 2
fi

mkdir -p "$OUTROOT"
"$PYTHON_BIN" verify_mpz_v9_14_install.py .
"$PYTHON_BIN" -m pytest -q \
  tests/test_event_remesh_czm_v914.py \
  tests/test_event_equilibrium_v914.py \
  tests/test_event_remesh_audit_v914.py \
  tests/test_mpz_v9_14_runner.py

EXTRA=()
[[ "$STOCHASTIC_EMISSION" == "1" ]] && EXTRA+=(--stochastic-emission) || EXTRA+=(--no-stochastic-emission)
[[ "$SKIP_EXISTING" == "1" ]] && EXTRA+=(--skip-existing) || EXTRA+=(--no-skip-existing)

echo "v9.14 event-remeshed full FEM/CZM gate"
echo "  classes=$CLASSES"
echo "  event statistics=$EVENT_STATISTICS"
echo "  stochastic emission=$STOCHASTIC_EMISSION"
echo "  adaptive event coordinate=absolute integrated hazard action"
echo "  maximum accepted dB=$ADAPTIVE_EVENT_TARGET"
echo "  one physical Arrhenius event per equilibrium solve"
echo "  conservative forward-patch remesh after every accepted event"
echo "  same-time/same-load FEM, MPZ-profile and J re-equilibration required"

"$PYTHON_BIN" run_mpz_v9_14_event_remesh_gate.py \
  --parameter-root "$PARAMETER_ROOT" \
  --outroot "$OUTROOT" \
  --seeds "$SEEDS" \
  --classes "$CLASSES" \
  --T-K "$T_K" \
  --target-extension-um "$TARGET_EXT_UM" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --n-stagger "$N_STAGGER" \
  --print-every "$PRINT_EVERY" \
  --adaptive-event-target "$ADAPTIVE_EVENT_TARGET" \
  --da-phys-um "$DA_PHYS_UM" \
  --mpz-length-um "$MPZ_LENGTH_UM" \
  --mpz-n-bins "$MPZ_N_BINS" \
  --crystal-theta-deg "$CRYSTAL_THETA_DEG" \
  --save-snapshots "$SAVE_SNAPSHOTS" \
  --snapshot-cols "$SNAPSHOT_COLS" \
  --snapshot-by-extension-um "$SNAPSHOT_BY_EXT_UM" \
  --event-statistics "$EVENT_STATISTICS" \
  --event-remesh-target-h-m "$EVENT_REMESH_TARGET_H_M" \
  --event-remesh-patch-radius-um "$EVENT_REMESH_PATCH_RADIUS_UM" \
  --event-remesh-max-edge-splits "$EVENT_REMESH_MAX_EDGE_SPLITS" \
  --event-remesh-target-edge-factor "$EVENT_REMESH_TARGET_EDGE_FACTOR" \
  --event-remesh-min-quality "$EVENT_REMESH_MIN_QUALITY" \
  "${EXTRA[@]}" \
  2>&1 | tee "$OUTROOT/driver.log"
