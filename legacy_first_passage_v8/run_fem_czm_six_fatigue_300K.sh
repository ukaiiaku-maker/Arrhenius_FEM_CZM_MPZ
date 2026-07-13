#!/usr/bin/env bash
set -euo pipefail

# FEM/CZM analogue of the six PF fatigue K-cycle simulations.
# This wrapper uses the existing v8 cyclic hazard/FEM driver and switches only
# the crack-geometry backend to adaptive/edge-split CZM.
#
# Typical use from Arrhenius_FEM_CZM project root:
#   MODE=smoke bash run_fem_czm_six_fatigue_300K.sh
#   MODE=pilot bash run_fem_czm_six_fatigue_300K.sh
#   MODE=full  bash run_fem_czm_six_fatigue_300K.sh

PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
cd "$PROJECT_ROOT"

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
if [[ -z "${PYTHON_BIN:-}" ]]; then
  # Prefer the real interpreter path. Some noninteractive shells expose `conda`
  # differently, and the older command-substitution resolver could leave
  # PYTHON_BIN empty, causing a bare "-m compileall"/command-not-found error.
  for cand in     "/opt/homebrew/Caskroom/miniconda/base/envs/$CONDA_ENV/bin/python"     "$HOME/miniconda3/envs/$CONDA_ENV/bin/python"     "$HOME/anaconda3/envs/$CONDA_ENV/bin/python"     "/opt/anaconda3/envs/$CONDA_ENV/bin/python"     "/opt/miniconda3/envs/$CONDA_ENV/bin/python"
  do
    if [[ -x "$cand" ]]; then
      PYTHON_BIN="$cand"
      break
    fi
  done
fi
if [[ -z "${PYTHON_BIN:-}" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3)"
  else
    PYTHON_BIN="$(command -v python)"
  fi
fi
if [[ -z "${PYTHON_BIN:-}" || ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: PYTHON_BIN is empty or not executable. Set PYTHON_BIN=/path/to/env/bin/python" >&2
  exit 2
fi

RUNNER="${RUNNER:-run_v8_compare_1d_2d_K_sweep.py}"
CASES_CSV="${CASES_CSV:-fem_czm_fatigue_cases.csv}"
OUTROOT="${OUTROOT:-runs/fem_czm_six_fatigue_300K}"
MODE="${MODE:-smoke}"   # smoke | pilot | full
CASE_FILTER="${CASE_FILTER:-}"
KLIST_OVERRIDE="${KLIST_OVERRIDE:-}"

T="${T:-300}"
R="${R:-0.1}"
FREQ="${FREQ:-1000}"
CYCLES_MAX="${CYCLES_MAX:-2e14}"
BLOCKS="${BLOCKS:-5000}"
BLOCK_CYCLES="${BLOCK_CYCLES:-1e5}"
MAX_BLOCK_CYCLES="${MAX_BLOCK_CYCLES:-inf}"
CYCLE_PHASES="${CYCLE_PHASES:-8}"

# FEM/CZM mesh/backend controls. Keep da_phys=5e-6 to match prior PF local points.
NX="${NX:-24}"
NY="${NY:-48}"
TIP_H_FINE="${TIP_H_FINE:-1e-6}"
TIP_RATIO="${TIP_RATIO:-1.25}"
DA_PHYS="${DA_PHYS:-5e-6}"
TARGET_EXT_UM="${TARGET_EXT_UM:-250}"
TARGET_DA_PER_BLOCK_UM="${TARGET_DA_PER_BLOCK_UM:-5}"
MAX_DA_PER_BLOCK_UM="${MAX_DA_PER_BLOCK_UM:-10}"
SNAPSHOTS="${SNAPSHOTS:-0}"
SNAPSHOT_COLS="${SNAPSHOT_COLS:-6}"
SNAPSHOT_BY_EXT_UM="${SNAPSHOT_BY_EXT_UM:-25}"
MAKE_2D_PLOTS="${MAKE_2D_PLOTS:-0}"

# Adaptive cycle-step controls used in the prior PF fatigue atlas.
TARGET_DB="${TARGET_DB:-0.01}"
TARGET_DN_STORE="${TARGET_DN_STORE:-0.01}"
TARGET_DN_EMIT="${TARGET_DN_EMIT:-0.20}"
TARGET_DN_MOBILE="${TARGET_DN_MOBILE:-0.20}"

CRACK_BACKEND="${CRACK_BACKEND:-adaptive_czm}"
CZM_MAX_ANGLE_ERROR_DEG="${CZM_MAX_ANGLE_ERROR_DEG:-35}"
# The current cyclic driver exposes only --czm-max-angle-error-deg among
# the adaptive-CZM geometry controls. Newer monotonic-driver-only flags
# are intentionally not passed here.

# To keep the first FEM fatigue reproduction close to the flat PF fatigue runs,
# leave branching off unless you explicitly pass branch options in EXTRA_ARGS.
EXTRA_ARGS="${EXTRA_ARGS:-}"

case "$MODE" in
  smoke|pilot|full) ;;
  *) echo "ERROR: MODE must be smoke, pilot, or full; got '$MODE'" >&2; exit 2 ;;
esac

echo "python: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PYCHECK'
import sys
print("Python", sys.version.split()[0])
PYCHECK

if [[ ! -f "$RUNNER" ]]; then
  echo "ERROR: could not find $RUNNER in $PWD" >&2
  echo "The FEM/CZM fatigue wrapper reuses the existing cyclic v8 driver." >&2
  echo "Copy run_v8_compare_1d_2d_K_sweep.py into this project root or set RUNNER=/path/to/it." >&2
  exit 2
fi
if [[ ! -f "$CASES_CSV" ]]; then
  echo "ERROR: could not find $CASES_CSV in $PWD" >&2
  exit 2
fi

# Check that the runner is the FEM/CZM-capable version. The architecture patch
# added --crack-backend to this driver. Refuse to run legacy PF geometry by accident.
if ! grep -q -- "--crack-backend" "$RUNNER"; then
  echo "ERROR: $RUNNER does not appear to support --crack-backend." >&2
  echo "Use the FEM/CZM-patched driver or apply the adaptive-CZM backend patch first." >&2
  exit 2
fi

"$PYTHON_BIN" -m compileall -q arrhenius_fracture
mkdir -p "$OUTROOT"
cp "$CASES_CSV" "$OUTROOT/fem_czm_fatigue_cases.csv"

plot_flag="--make-2d-plots"
if [[ "$MAKE_2D_PLOTS" == "0" || "$MAKE_2D_PLOTS" == "false" ]]; then
  plot_flag="--no-make-2d-plots"
fi

should_run_case() {
  local label="$1"
  if [[ -z "$CASE_FILTER" ]]; then
    if [[ "$MODE" == "smoke" ]]; then
      case "$label" in
        FCC_like_case29|plastic_shielded_case64_M1|steep_cleavage_case35) return 0 ;;
        *) return 1 ;;
      esac
    fi
    return 0
  fi
  local normalized=" ${CASE_FILTER//,/ } "
  [[ "$normalized" == *" $label "* ]]
}

run_one_case() {
  local label="$1"
  local source_case="$2"
  local material_class="$3"
  local full_klist="$4"
  local smoke_klist="$5"
  local pilot_klist="$6"
  local G00="$7"
  local sigc="$8"
  local expa="$9"
  local expn="${10}"
  local floor="${11}"
  local emitS="${12}"
  local peierlsS="${13}"
  local taylorS="${14}"

  if ! should_run_case "$label"; then
    echo "=== Skipping $label ==="
    return 0
  fi

  local klist
  if [[ -n "$KLIST_OVERRIDE" ]]; then
    klist="$KLIST_OVERRIDE"
  elif [[ "$MODE" == "full" ]]; then
    klist="$full_klist"
  elif [[ "$MODE" == "pilot" ]]; then
    klist="$pilot_klist"
  else
    klist="$smoke_klist"
  fi

  echo
  echo "======================================================================"
  echo "FEM/CZM fatigue case: $label"
  echo "class: $material_class"
  echo "mode=$MODE; Kmax=[$klist]; T=$T K; R=$R; f=$FREQ Hz"
  echo "target_ext=$TARGET_EXT_UM um; cycles_max=$CYCLES_MAX; blocks=$BLOCKS"
  echo "backend=$CRACK_BACKEND"
  echo "======================================================================"

  mkdir -p "$OUTROOT/$label"

  # shellcheck disable=SC2086
  "$PYTHON_BIN" "$RUNNER" \
    --out "$OUTROOT/$label" \
    --Kmax-MPa-sqrt-m $klist \
    --T "$T" \
    --R "$R" \
    --frequency-Hz "$FREQ" \
    --blocks "$BLOCKS" \
    --cycles-max "$CYCLES_MAX" \
    --block-cycles "$BLOCK_CYCLES" \
    --max-block-cycles "$MAX_BLOCK_CYCLES" \
    --cycle-block-mode hazard_limited \
    --target-dB "$TARGET_DB" \
    --target-dN-store "$TARGET_DN_STORE" \
    --target-dN-emit "$TARGET_DN_EMIT" \
    --target-dN-mobile "$TARGET_DN_MOBILE" \
    --storage-model escape_limited \
    --calibrate-2d-K \
    --K-calib-iters 3 \
    --K-calib-tol 5e-3 \
    --no-stop-after-first-2d-fire \
    --cyclic-mechanics-phases "$CYCLE_PHASES" \
    --nx "$NX" --ny "$NY" \
    --tip-h-fine "$TIP_H_FINE" \
    --tip-ratio "$TIP_RATIO" \
    --da-phys "$DA_PHYS" \
    --target-da-per-block-um "$TARGET_DA_PER_BLOCK_UM" \
    --target-crack-extension-um "$TARGET_EXT_UM" \
    --snapshot-by-crack-extension-um "$SNAPSHOT_BY_EXT_UM" \
    --max-da-per-block-um "$MAX_DA_PER_BLOCK_UM" \
    --save-snapshots "$SNAPSHOTS" \
    --snapshot-cols "$SNAPSHOT_COLS" \
    $plot_flag \
    --min-global-forward 0.05 \
    --crack-backend "$CRACK_BACKEND" \
    --czm-max-angle-error-deg "$CZM_MAX_ANGLE_ERROR_DEG" \
    --cleave-barrier-kind exp_floor \
    --cleave-exp-T-mode mu_scale \
    --cleave-G00-eV "$G00" \
    --cleave-sigc0-GPa "$sigc" \
    --cleave-exp-a "$expa" \
    --cleave-exp-n "$expn" \
    --cleave-floor-frac "$floor" \
    --emit-energy-scale 0.75 \
    --emit-entropy-scale "$emitS" \
    --peierls-energy-scale 0.00375 \
    --peierls-entropy-scale "$peierlsS" \
    --peierls-stress-scale 1.0 \
    --taylor-energy-scale 0.015 \
    --taylor-entropy-scale "$taylorS" \
    --taylor-stress-scale 1.0 \
    $EXTRA_ARGS
}

# Use Python's CSV parser so quoted K lists are handled correctly under bash 3.2.
while IFS=$'\t' read -r label source_case material_class full_klist smoke_klist pilot_klist G00 sigc expa expn floor emitS peierlsS taylorS; do
  run_one_case "$label" "$source_case" "$material_class" "$full_klist" "$smoke_klist" "$pilot_klist" "$G00" "$sigc" "$expa" "$expn" "$floor" "$emitS" "$peierlsS" "$taylorS"
done < <("$PYTHON_BIN" - <<'PY' "$CASES_CSV"
import csv, sys
with open(sys.argv[1], newline='') as f:
    r = csv.DictReader(f)
    for row in r:
        fields = [
            row['case_label'], row['source_case'], row['material_response_class'],
            row['full_klist'], row['smoke_klist'], row['pilot_klist'],
            row['cleave_G00_eV'], row['cleave_sigc0_GPa'], row['cleave_exp_a'],
            row['cleave_exp_n'], row['cleave_floor_frac'], row['emit_entropy_scale'],
            row['peierls_entropy_scale'], row['taylor_entropy_scale'],
        ]
        print('\t'.join(fields))
PY
)

if [[ "${RUN_ANALYSIS:-1}" != "0" && -f analyze_fem_czm_fatigue_outputs.py ]]; then
  "$PYTHON_BIN" analyze_fem_czm_fatigue_outputs.py --root "$OUTROOT" --R "$R" --cycles-max "$CYCLES_MAX" --target-crack-extension-um "$TARGET_EXT_UM"
fi

echo
echo "DONE. Output root: $OUTROOT"
