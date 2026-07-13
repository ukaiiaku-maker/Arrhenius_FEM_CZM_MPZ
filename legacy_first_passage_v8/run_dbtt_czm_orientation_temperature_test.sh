#!/usr/bin/env bash
set -euo pipefail

THETA=${THETA:-30}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000"}
OUTROOT=${OUTROOT:-runs/dbtt_czm_theta${THETA}_temperature_test}
BACKEND=${CRACK_BACKEND:-adaptive_czm}
MAX_JOBS=${MAX_JOBS:-2}
FORCE=${FORCE:-0}
CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-0}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-4}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-0}
LONG_GROWTH=${LONG_GROWTH:-0}
TARGET_EXT_UM=${TARGET_EXT_UM:-1200}
STEPS=${STEPS:-4500}
LONG_STEPS=${LONG_STEPS:-12000}
MAX_FRONTS=${MAX_FRONTS:-32}
RETIRE_STAGNANT_BRANCHES=${RETIRE_STAGNANT_BRANCHES:-0}
COALESCE_CRACKS=${COALESCE_CRACKS:-1}

# Canonical fracture-response regime controls.  Defaults reproduce the DBTT
# class; outer sweep drivers can override these without duplicating the full
# mechanics/hazard command.
CLEAVE_H0_EV=${CLEAVE_H0_EV:-6.0}
CLEAVE_SHIELD_CHI=${CLEAVE_SHIELD_CHI:-0.60}
N_SAT=${N_SAT:-2000}
FRACTURE_CLASS=${FRACTURE_CLASS:-DBTT}

# Optional branch-birth controls.  PERMISSIVE_BRANCHING=1 enables the branch
# clock and applies a deliberately permissive, but still mechanically gated,
# preset.  Individual values may be overridden from the environment.
PERMISSIVE_BRANCHING=${PERMISSIVE_BRANCHING:-0}
ENABLE_BRANCHING=${ENABLE_BRANCHING:-0}
if [[ "$PERMISSIVE_BRANCHING" == "1" ]]; then
  ENABLE_BRANCHING=1
  BRANCH_FP_MIN_RATIO=${BRANCH_FP_MIN_RATIO:-0.75}
  BRANCH_CLOCK_TARGET=${BRANCH_CLOCK_TARGET:-0.50}
  BRANCH_SECONDARY_MIN_K_RATIO=${BRANCH_SECONDARY_MIN_K_RATIO:-0.75}
  BRANCH_SPACING=${BRANCH_SPACING:-5.0}
else
  BRANCH_FP_MIN_RATIO=${BRANCH_FP_MIN_RATIO:-0.95}
  BRANCH_CLOCK_TARGET=${BRANCH_CLOCK_TARGET:-1.0}
  BRANCH_SECONDARY_MIN_K_RATIO=${BRANCH_SECONDARY_MIN_K_RATIO:-0.85}
  BRANCH_SPACING=${BRANCH_SPACING:-10.0}
fi

# Interpreter precedence:
#   1. explicit PYTHON_BIN
#   2. Python from CONDA_ENV (default: arrhenius-fem-czm)
# This deliberately avoids silently using a broken or mutable base environment.
if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_EXE="$PYTHON_BIN"
elif command -v conda >/dev/null 2>&1 && conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV"; then
  PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 \
    | tr -d '\r' \
    | awk 'NF {last=$0} END {print last}')"
else
  echo "ERROR: isolated Conda environment '$CONDA_ENV' was not found." >&2
  echo "Create it once with:" >&2
  echo "  bash setup_fem_czm_env.sh" >&2
  echo "or override explicitly with PYTHON_BIN=/path/to/python." >&2
  exit 2
fi

if [[ -z "$PYTHON_EXE" || ! -x "$PYTHON_EXE" ]]; then
  echo "ERROR: no usable Python interpreter found." >&2
  exit 2
fi

mkdir -p "$OUTROOT"

echo "=== DBTT adaptive-CZM orientation sweep ==="
echo "project:  $PWD"
echo "python:   $PYTHON_EXE"
"$PYTHON_EXE" --version
echo "theta:    ${THETA} deg"
echo "temps:    ${TEMPS}"
echo "backend:  ${BACKEND}"
echo "conda_env:${CONDA_ENV}"
echo "max_jobs: ${MAX_JOBS}"
echo "outroot:  ${OUTROOT}"
echo "snapshots: ${SAVE_SNAPSHOTS} (cols=${SNAPSHOT_COLS}, by_ext_um=${SNAPSHOT_BY_EXT_UM})"
echo "long_growth: ${LONG_GROWTH} (target_ext_um=${TARGET_EXT_UM}, steps=${STEPS}, long_steps=${LONG_STEPS})"
echo "coalescence: ${COALESCE_CRACKS}"
echo "fracture_class: ${FRACTURE_CLASS} (H0=${CLEAVE_H0_EV} eV, chi=${CLEAVE_SHIELD_CHI}, n_sat=${N_SAT})"
if [[ "$ENABLE_BRANCHING" == "1" ]]; then
  echo "branching: ON (permissive=${PERMISSIVE_BRANCHING}, fp_min_ratio=${BRANCH_FP_MIN_RATIO}, clock_target=${BRANCH_CLOCK_TARGET}, secondary_K_ratio=${BRANCH_SECONDARY_MIN_K_RATIO}, spacing_da=${BRANCH_SPACING}, max_fronts=${MAX_FRONTS}, retire_stagnant=${RETIRE_STAGNANT_BRANCHES})"
else
  echo "branching: OFF"
fi

# Fail before creating a silent batch of broken jobs.
"$PYTHON_EXE" - <<'PY'
import sys
import numpy
import scipy
import arrhenius_fracture
from arrhenius_fracture.config import FractureBarrier
print(f"preflight OK: Python {sys.version.split()[0]}, numpy {numpy.__version__}, scipy {scipy.__version__}")
PY

run_one() {
  local T="$1"
  local OUT="$OUTROOT/T${T}_th${THETA}"
  local LOG="$OUT/run.log"

  local SNAPSHOT_PNG="$OUT/field_snapshots_${T}K.png"

  if [[ -f "$OUT/summary.json" && "$FORCE" != "1" ]]; then
    if [[ "$LONG_GROWTH" == "1" && ! -f "$OUT/.long_growth_complete" ]]; then
      echo "RERUN T=${T} K: summary exists but long-growth target is not marked complete"
    elif [[ "$SAVE_SNAPSHOTS" -gt 0 && ! -f "$SNAPSHOT_PNG" ]]; then
      echo "RERUN T=${T} K: summary exists but requested snapshot PNG is missing"
    else
      echo "SKIP T=${T} K: existing completed outputs"
      return 0
    fi
  fi

  mkdir -p "$OUT"
  echo "START T=${T} K -> $OUT"

  local plot_args=()
  local growth_args=()
  local steps_run="$STEPS"
  local snapshot_ext_run="$SNAPSHOT_BY_EXT_UM"
  local branch_args=()
  local coalescence_args=()
  if [[ "$COALESCE_CRACKS" == "1" ]]; then
    coalescence_args+=(--coalesce-cracks)
  else
    coalescence_args+=(--no-coalesce-cracks)
  fi
  if [[ "$ENABLE_BRANCHING" == "1" ]]; then
    branch_args+=(
      --crystal-branch
      --branch-fp-min-ratio "$BRANCH_FP_MIN_RATIO"
      --branch-clock-target "$BRANCH_CLOCK_TARGET"
      --branch-secondary-min-K-ratio "$BRANCH_SECONDARY_MIN_K_RATIO"
      --branch-spacing "$BRANCH_SPACING"
    )
    if [[ "$RETIRE_STAGNANT_BRANCHES" == "1" ]]; then
      branch_args+=(--retire-stagnant-branches)
    fi
  fi

  if [[ "$LONG_GROWTH" == "1" ]]; then
    steps_run="$LONG_STEPS"
    growth_args+=(--target-crack-extension-um "$TARGET_EXT_UM")
    # Long-growth output is most informative when snapshots are tied to physical
    # crack extension rather than nominal load step.  Respect an explicit user
    # value; otherwise use 100 um spacing.
    if [[ "$snapshot_ext_run" == "0" || "$snapshot_ext_run" == "0.0" ]]; then
      snapshot_ext_run=100
    fi
  else
    growth_args+=(--stop-after-first-fire)
  fi

  if [[ "$SAVE_SNAPSHOTS" -gt 0 ]]; then
    plot_args+=(--save-snapshots "$SAVE_SNAPSHOTS" --snapshot-cols "$SNAPSHOT_COLS")
    if [[ "$snapshot_ext_run" != "0" && "$snapshot_ext_run" != "0.0" ]]; then
      plot_args+=(--snapshot-by-crack-extension-um "$snapshot_ext_run")
    fi
  else
    plot_args+=(--save-snapshots 0 --no-plots)
  fi

  if "$PYTHON_EXE" -m arrhenius_fracture.sharp_front \
      --mode 2d --nx 12 --ny 24 \
      --tip-h-fine 5e-6 --tip-ratio 1.30 \
      --dU 2e-7 --dt 8.4 --steps "$steps_run" --n-stagger 2 \
      "${plot_args[@]}" "${growth_args[@]}" "${branch_args[@]}" "${coalescence_args[@]}" --print-every 500 \
      --crystal-aniso --crystal-compete --crystal-material branchy \
      --cleave-gamma-aniso 2.0 \
      --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
      --multihit-m 3 --multihit-tau 1e-6 \
      --emb-sat-frac 1 \
      --adaptive-events --adaptive-event-target 0.35 --adaptive-min-frac 1e-8 --adaptive-grow 4.0 \
      --max-fronts "$MAX_FRONTS" --da-phys 5e-6 \
      --j-decomposition cluster --rJ-cluster 20e-6 --rJ-outer 25e-6 \
      --cleave-H0-eV "$CLEAVE_H0_EV" --cleave-shield-chi "$CLEAVE_SHIELD_CHI" --n-sat "$N_SAT" \
      --temperatures "$T" --crystal-theta-deg "$THETA" \
      --crack-backend "$BACKEND" --czm-max-angle-error-deg 35 \
      --out "$OUT" > "$LOG" 2>&1; then
    echo "DONE  T=${T} K"
    if [[ "$LONG_GROWTH" == "1" ]]; then
      if grep -Eq "reached target crack extension|ligament severed" "$LOG"; then
        touch "$OUT/.long_growth_complete"
      else
        rm -f "$OUT/.long_growth_complete"
        echo "WARNING T=${T} K completed solver steps without reaching the requested long-growth target" >&2
      fi
    fi
  else
    local rc=$?
    echo "FAILED T=${T} K (exit ${rc})" >&2
    echo "----- tail: $LOG -----" >&2
    tail -n 60 "$LOG" >&2 || true
    echo "----- end tail -----" >&2
    return "$rc"
  fi
}

export -f run_one
export THETA OUTROOT BACKEND FORCE PYTHON_EXE CONDA_ENV SAVE_SNAPSHOTS SNAPSHOT_COLS SNAPSHOT_BY_EXT_UM LONG_GROWTH TARGET_EXT_UM STEPS LONG_STEPS PERMISSIVE_BRANCHING ENABLE_BRANCHING BRANCH_FP_MIN_RATIO BRANCH_CLOCK_TARGET BRANCH_SECONDARY_MIN_K_RATIO BRANCH_SPACING MAX_FRONTS RETIRE_STAGNANT_BRANCHES COALESCE_CRACKS CLEAVE_H0_EV CLEAVE_SHIELD_CHI N_SAT FRACTURE_CLASS

# `bash -c`, not `bash -lc`: a login shell can reset PATH and bypass Conda.
if ! printf '%s\n' $TEMPS | xargs -n1 -P "$MAX_JOBS" bash -c 'run_one "$1"' _; then
  echo "ERROR: one or more temperature jobs failed. See the FAILED block above and per-case run.log files." >&2
  exit 1
fi

"$PYTHON_EXE" - <<'PY'
import json, csv, glob, os, math
import numpy as np

root=os.environ['OUTROOT']
rows=[]
for fn in glob.glob(os.path.join(root,'T*_th*','summary.json')):
    with open(fn) as f:
        data=json.load(f)
    if not data:
        continue
    d=data[0]
    T=d['T']
    path=glob.glob(os.path.join(os.path.dirname(fn),f'crack_path_{int(T)}K.csv'))
    inserted_angle=float('nan')
    if path:
        a=np.loadtxt(path[0],delimiter=',',skiprows=1)
        if a.ndim==2 and len(a)>=2:
            dx,dy=a[-1]-a[-2]
            inserted_angle=math.degrees(math.atan2(dy,dx))
    bdiag=glob.glob(os.path.join(os.path.dirname(fn),f'branch_diagnostics_{int(T):04d}K.csv'))
    hazard_angle=float('nan')
    if bdiag:
        with open(bdiag[0], newline='') as bf:
            rr=list(csv.DictReader(bf))
            if rr:
                hazard_angle=float(rr[-1]['angle1_deg'])
    rows.append({
        'T_K':T,
        'Kc_first_MPa_sqrt_m':d['Kc_first_MPa_sqrt_m'],
        'N_em_final':d['N_em_final'],
        'deflection_deg_summary':d['deflection_deg'],
        'hazard_selected_angle_deg':hazard_angle,
        'inserted_segment_angle_deg':inserted_angle,
        'angle_error_deg':abs(inserted_angle-hazard_angle),
        'mode':d['mode'],
    })

rows.sort(key=lambda r:r['T_K'])
if not rows:
    raise SystemExit(f"ERROR: no completed summary.json files found under {root}")

out=os.path.join(root,'dbtt_czm_temperature_summary.csv')
with open(out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"WROTE {out}")
for r in rows:
    print(r)
PY
