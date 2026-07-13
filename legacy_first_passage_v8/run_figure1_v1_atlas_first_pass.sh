#!/usr/bin/env bash
set -euo pipefail

# Figure 1 V1 atlas, first pass: Panels A, C, and F.
#
# This workflow does not alter solver physics. It calls the existing V1 source
# workflows, preserves their raw CSV outputs, then assembles publication-style
# provisional panels plus panel-specific replot-ready CSVs.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
OUTROOT="${OUTROOT:-runs/figure1_v1_atlas_first_pass}"

FORWARD_RUNNER="${FORWARD_RUNNER:-run_forward_barrier_KT_prediction.py}"
REFINED_DRIVER="${REFINED_DRIVER:-run_refined_two_barrier_dbtt_fatigue_strength.sh}"
BUILDER="${BUILDER:-$HERE/build_figure1_v1_atlas_first_pass.py}"
MANIFEST="${MANIFEST:-$HERE/figure1_case_manifest.csv}"
CASE_TABLE="${CASE_TABLE:-selected_v1_temperature_cases_corrected.csv}"
MAP_RUNNER="${MAP_RUNNER:-run_v1_two_barrier_dbtt_fatigue_map_fixed.py}"

RUN_PANEL_A_SOURCE="${RUN_PANEL_A_SOURCE:-1}"
RUN_CORE_SOURCE="${RUN_CORE_SOURCE:-1}"

PRIMARY_RATE="${PRIMARY_RATE:-1e-10}"
SENSITIVITY_RATE="${SENSITIVITY_RATE:-1e-12}"
S_EMIT="${S_EMIT:--40}"
S_CLEAVE="${S_CLEAVE:-0}"
CORE_TEMPS="${CORE_TEMPS:-100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000}"

PANEL_A_DIR="$OUTROOT/source_panel_A"
CORE_DIR="$OUTROOT/source_core"
ASSEMBLY_DIR="$OUTROOT/assembled_first_pass"

mkdir -p "$OUTROOT"

# Preflight the active environment before expensive work.
"$PYTHON" - <<'PY'
import sys
import numpy, pandas, scipy, matplotlib
print("python:", sys.executable)
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("scipy:", scipy.__version__)
print("matplotlib:", matplotlib.__version__)
PY

if [[ ! -f "$BUILDER" ]]; then
  echo "ERROR: builder not found: $BUILDER" >&2
  exit 2
fi
if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2
  exit 2
fi

if [[ "$RUN_PANEL_A_SOURCE" == "1" ]]; then
  if [[ ! -f "$FORWARD_RUNNER" ]]; then
    echo "ERROR: Panel A runner not found: $FORWARD_RUNNER" >&2
    exit 2
  fi
  mkdir -p "$PANEL_A_DIR"
  echo "=== Panel A source: four canonical Kc(T) classes ==="
  "$PYTHON" "$FORWARD_RUNNER" \
    --out "$PANEL_A_DIR" \
    --temperatures 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000 1025 1050 1075 1100 1125 1150 1175 1200 1225 1250 1275 1300 \
    --regimes ceramic peak weakT dbtt \
    --Kdot "${PANEL_A_KDOT:-0.02}" \
    --Kmax "${PANEL_A_KMAX:-80}" \
    --dt "${PANEL_A_DT:-1.0}"
fi

if [[ "$RUN_CORE_SOURCE" == "1" ]]; then
  if [[ ! -f "$REFINED_DRIVER" ]]; then
    echo "ERROR: refined V1 driver not found: $REFINED_DRIVER" >&2
    exit 2
  fi
  echo "=== Core source: six-case V1 Kc(T) + adaptive fatigue thresholds ==="
  OUTROOT="$CORE_DIR" \
  CASE_TABLE="$CASE_TABLE" \
  MAP_RUNNER="$MAP_RUNNER" \
  RUN_EXISTING_ANALYSIS=0 \
  RUN_CORE_ADAPTIVE=1 \
  RUN_EMISSION_STUDY=0 \
  RUN_STRENGTH=0 \
  RUN_CORRELATION=0 \
  PRIMARY_RATE="$PRIMARY_RATE" \
  SENSITIVITY_RATE="$SENSITIVITY_RATE" \
  CORE_TEMPS="$CORE_TEMPS" \
  bash "$REFINED_DRIVER"
fi

# Resolve Panel A CSV name across the two currently used forward-runner variants.
PANEL_A_CSV=""
for p in \
  "$PANEL_A_DIR/forward_barrier_KT_prediction.csv" \
  "$PANEL_A_DIR/input_parameter_KT_prediction.csv"
do
  if [[ -f "$p" ]]; then PANEL_A_CSV="$p"; break; fi
done
if [[ -z "$PANEL_A_CSV" ]]; then
  echo "ERROR: could not find Panel A prediction CSV under $PANEL_A_DIR" >&2
  exit 2
fi

FATIGUE_CSV="$CORE_DIR/core_adaptive/fatigue_adaptive_points.csv"
MONO_CSV="$CORE_DIR/core_adaptive/monotonic_adaptive_points.csv"
THRESH_CSV="$CORE_DIR/core_adaptive_analysis/rate_defined_thresholds.csv"

for p in "$FATIGUE_CSV" "$MONO_CSV" "$THRESH_CSV"; do
  if [[ ! -f "$p" ]]; then
    echo "ERROR: required core source missing: $p" >&2
    exit 2
  fi
done

echo "=== Assembling Panels A, C, F and replot-ready CSVs ==="
"$PYTHON" "$BUILDER" \
  --panel-a-csv "$PANEL_A_CSV" \
  --fatigue-csv "$FATIGUE_CSV" \
  --monotonic-csv "$MONO_CSV" \
  --threshold-csv "$THRESH_CSV" \
  --manifest "$MANIFEST" \
  --out "$ASSEMBLY_DIR" \
  --S-emit "$S_EMIT" \
  --S-cleave "$S_CLEAVE" \
  --panel-c-T "${PANEL_C_T:-300}" \
  --metric-T-ref "${METRIC_T_REF:-300}" \
  --metric-T-hi-ref "${METRIC_T_HI_REF:-900}" \
  --primary-rate "$PRIMARY_RATE"

echo
echo "Figure 1 first-pass atlas complete."
echo "Raw Panel A source: $PANEL_A_DIR"
echo "Raw six-case source: $CORE_DIR"
echo "Assembled panels/data: $ASSEMBLY_DIR"
