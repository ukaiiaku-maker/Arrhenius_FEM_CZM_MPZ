#!/usr/bin/env bash
set -euo pipefail

# Production long-growth sweep over the four canonical fracture-temperature
# response classes used by the previous PF/regime calculations.
#
# Defaults:
#   classes: ceramic peak weakT DBTT
#   T:       300..1200 K in 100 K increments
#   crystal: 45 deg
#   morphology: PF-like branching + exact coalescence
#
# Useful overrides:
#   THETA=30
#   CLASSES="ceramic DBTT"
#   TEMPS="300 400 500"
#   TARGET_EXT_UM=750
#   MAX_JOBS=1
#   FORCE=1
#   OUTROOT=runs/four_class_fracture_theta45

THETA=${THETA:-45}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100 1200"}
CLASSES=${CLASSES:-"ceramic peak weakT DBTT"}
OUTROOT=${OUTROOT:-runs/four_class_fracture_theta${THETA}_pf_like}
MAX_JOBS=${MAX_JOBS:-1}
FORCE=${FORCE:-0}
TARGET_EXT_UM=${TARGET_EXT_UM:-750}
LONG_STEPS=${LONG_STEPS:-20000}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-12}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-6}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-75}
CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}

# Keep morphology controls common across all four classes.
BRANCH_FP_MIN_RATIO=${BRANCH_FP_MIN_RATIO:-0.87}
BRANCH_CLOCK_TARGET=${BRANCH_CLOCK_TARGET:-0.90}
BRANCH_SECONDARY_MIN_K_RATIO=${BRANCH_SECONDARY_MIN_K_RATIO:-0.85}
BRANCH_SPACING=${BRANCH_SPACING:-10.0}
MAX_FRONTS=${MAX_FRONTS:-3}
RETIRE_STAGNANT_BRANCHES=${RETIRE_STAGNANT_BRANCHES:-1}
COALESCE_CRACKS=${COALESCE_CRACKS:-1}

mkdir -p "$OUTROOT"

regime_values() {
  case "$1" in
    ceramic) echo "2.6 0.00 inf" ;;
    peak)    echo "3.6 0.10 inf" ;;
    weakT)   echo "4.0 0.20 1500" ;;
    DBTT)    echo "6.0 0.60 2000" ;;
    *)
      echo "ERROR: unknown fracture class '$1'" >&2
      return 2
      ;;
  esac
}

cat > "$OUTROOT/four_class_sweep_config.txt" <<EOF
Four-class adaptive-CZM fracture temperature sweep
THETA=${THETA}
TEMPS=${TEMPS}
CLASSES=${CLASSES}
TARGET_EXT_UM=${TARGET_EXT_UM}
LONG_STEPS=${LONG_STEPS}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM}
MAX_FRONTS=${MAX_FRONTS}
COALESCE_CRACKS=${COALESCE_CRACKS}

Canonical class parameters:
  ceramic : H0=2.6 eV, chi=0.00, n_sat=inf
  peak    : H0=3.6 eV, chi=0.10, n_sat=inf
  weakT   : H0=4.0 eV, chi=0.20, n_sat=1500
  DBTT    : H0=6.0 eV, chi=0.60, n_sat=2000

Shared PF-like morphology controls:
  branch_fp_min_ratio=${BRANCH_FP_MIN_RATIO}
  branch_clock_target=${BRANCH_CLOCK_TARGET}
  branch_secondary_min_K_ratio=${BRANCH_SECONDARY_MIN_K_RATIO}
  branch_spacing_da=${BRANCH_SPACING}
  max_active_fronts=${MAX_FRONTS}
  retire_stagnant_branches=${RETIRE_STAGNANT_BRANCHES}
EOF

echo "=== Four-class fracture sweep ==="
echo "theta:     ${THETA} deg"
echo "classes:   ${CLASSES}"
echo "temps:     ${TEMPS}"
echo "target:    ${TARGET_EXT_UM} um"
echo "outroot:   ${OUTROOT}"
echo "max_jobs:  ${MAX_JOBS} temperature jobs per class"

FAILED_CLASSES=()
for CLASS in $CLASSES; do
  read -r H0 CHI NSAT <<< "$(regime_values "$CLASS")"
  CLASS_OUT="$OUTROOT/$CLASS"
  echo
  echo "=== CLASS ${CLASS}: H0=${H0} eV, chi=${CHI}, n_sat=${NSAT} ==="

  if ! THETA="$THETA" \
  TEMPS="$TEMPS" \
  OUTROOT="$CLASS_OUT" \
  MAX_JOBS="$MAX_JOBS" \
  FORCE="$FORCE" \
  TARGET_EXT_UM="$TARGET_EXT_UM" \
  LONG_STEPS="$LONG_STEPS" \
  SAVE_SNAPSHOTS="$SAVE_SNAPSHOTS" \
  SNAPSHOT_COLS="$SNAPSHOT_COLS" \
  SNAPSHOT_BY_EXT_UM="$SNAPSHOT_BY_EXT_UM" \
  CONDA_ENV="$CONDA_ENV" \
  BRANCH_FP_MIN_RATIO="$BRANCH_FP_MIN_RATIO" \
  BRANCH_CLOCK_TARGET="$BRANCH_CLOCK_TARGET" \
  BRANCH_SECONDARY_MIN_K_RATIO="$BRANCH_SECONDARY_MIN_K_RATIO" \
  BRANCH_SPACING="$BRANCH_SPACING" \
  MAX_FRONTS="$MAX_FRONTS" \
  RETIRE_STAGNANT_BRANCHES="$RETIRE_STAGNANT_BRANCHES" \
  COALESCE_CRACKS="$COALESCE_CRACKS" \
  FRACTURE_CLASS="$CLASS" \
  CLEAVE_H0_EV="$H0" \
  CLEAVE_SHIELD_CHI="$CHI" \
  N_SAT="$NSAT" \
  bash run_dbtt_czm_pf_like_branching.sh; then
    echo "WARNING: class ${CLASS} had one or more failed temperature cases; continuing with remaining classes" >&2
    FAILED_CLASSES+=("$CLASS")
  fi
done

# Build one compact cross-class summary from the per-case summary JSON files.
if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_EXE="$PYTHON_BIN"
else
  PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 \
    | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
fi

OUTROOT="$OUTROOT" CLASSES="$CLASSES" "$PYTHON_EXE" - <<'PY'
import csv, glob, json, os
from pathlib import Path

root = Path(os.environ['OUTROOT'])
classes = os.environ['CLASSES'].split()
rows = []
for cls in classes:
    for fn in glob.glob(str(root / cls / 'T*_th*' / 'summary.json')):
        p = Path(fn)
        with p.open() as f:
            data = json.load(f)
        if not data:
            continue
        d = data[0]
        rows.append({
            'class': cls,
            'T_K': float(d['T']),
            'Kc_first_MPa_sqrt_m': d.get('Kc_first_MPa_sqrt_m'),
            'N_em_final': d.get('N_em_final'),
            'deflection_deg': d.get('deflection_deg'),
            'n_active_fronts_final': d.get('n_active_fronts_final'),
            'n_coalesced': d.get('n_coalesced'),
            'mode': d.get('mode'),
            'case_dir': str(p.parent.relative_to(root)),
        })
rows.sort(key=lambda r: (classes.index(r['class']), r['T_K']))
out = root / 'four_class_temperature_summary.csv'
fields = ['class','T_K','Kc_first_MPa_sqrt_m','N_em_final','deflection_deg',
          'n_active_fronts_final','n_coalesced','mode','case_dir']
with out.open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader(); w.writerows(rows)
print(f'WROTE {out} ({len(rows)} completed class-temperature cases)')

# Reader-facing class comparison plots.
if rows:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.0, 4.8))
        for cls in classes:
            rr = [r for r in rows if r['class'] == cls]
            if rr:
                ax.plot([r['T_K'] for r in rr],
                        [r['Kc_first_MPa_sqrt_m'] for r in rr],
                        marker='o', label=cls)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel(r'First-passage $K_c$ (MPa$\sqrt{m}$)')
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(root / 'four_class_Kc_vs_T.png', dpi=220)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7.0, 4.8))
        for cls in classes:
            rr = [r for r in rows if r['class'] == cls]
            if rr:
                ax.plot([r['T_K'] for r in rr],
                        [r['N_em_final'] for r in rr],
                        marker='o', label=cls)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel(r'$N_{em}$ at end of run')
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(root / 'four_class_Nem_vs_T.png', dpi=220)
        plt.close(fig)
        print('WROTE four_class_Kc_vs_T.png and four_class_Nem_vs_T.png')
    except Exception as exc:
        print(f'WARNING: comparison plotting failed: {exc}')
PY

if (( ${#FAILED_CLASSES[@]} > 0 )); then
  printf '%s\n' "${FAILED_CLASSES[@]}" > "$OUTROOT/failed_classes.txt"
  echo "WARNING: failed classes: ${FAILED_CLASSES[*]}" >&2
  echo "Completed/partial results and combined summaries were retained." >&2
  exit 1
fi

echo "=== Four-class sweep complete ==="
