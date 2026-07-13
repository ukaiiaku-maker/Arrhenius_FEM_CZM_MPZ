#!/usr/bin/env bash
set -euo pipefail

# Refined two-barrier study:
#   0. Re-analyze an existing corrected map with rate-defined thresholds.
#   1. Adaptive Kc(T) / DeltaK_th(T) study for the six canonical response classes.
#   2. Emission-entropy sweep for the fatigue/strength-anomaly connection.
#   3. Arrhenius strength-anomaly calculation from the same emission barrier.
#   4. Correlation/residual analysis.
#
# This script does not overwrite the original corrected map.  All new output is
# written below OUTROOT.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
OUTROOT="${OUTROOT:-runs/refined_two_barrier_dbtt_fatigue_strength}"
EXISTING_MAP="${EXISTING_MAP:-runs/v1_two_barrier_dbtt_fatigue_map_corrected}"
CASE_TABLE="${CASE_TABLE:-selected_v1_temperature_cases_corrected.csv}"
MAP_RUNNER="${MAP_RUNNER:-}"

RUN_EXISTING_ANALYSIS="${RUN_EXISTING_ANALYSIS:-1}"
RUN_CORE_ADAPTIVE="${RUN_CORE_ADAPTIVE:-1}"
RUN_EMISSION_STUDY="${RUN_EMISSION_STUDY:-1}"
RUN_STRENGTH="${RUN_STRENGTH:-1}"
RUN_CORRELATION="${RUN_CORRELATION:-1}"

# Production defaults. Override from the environment as needed.
CYCLES_MAX="${CYCLES_MAX:-2e14}"
MAX_BLOCKS="${MAX_BLOCKS:-10000}"
PRIMARY_RATE="${PRIMARY_RATE:-1e-10}"
SENSITIVITY_RATE="${SENSITIVITY_RATE:-1e-12}"
CORE_TEMPS="${CORE_TEMPS:-250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000}"
EMISSION_TEMPS="${EMISSION_TEMPS:-250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000}"
STRENGTH_TEMPS="${STRENGTH_TEMPS:-250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000}"
EMISSION_ENTROPIES="${EMISSION_ENTROPIES:--20 -25 -30 -35 -40 -45 -50 -55 -60}"
STRENGTH_RATES="${STRENGTH_RATES:-1e-4 1e-2 1 100}"

mkdir -p "$OUTROOT"

if [[ -z "$MAP_RUNNER" ]]; then
  if [[ -f run_v1_two_barrier_dbtt_fatigue_map_fixed.py ]]; then
    MAP_RUNNER="run_v1_two_barrier_dbtt_fatigue_map_fixed.py"
  elif [[ -f run_v1_two_barrier_dbtt_fatigue_map.py ]]; then
    MAP_RUNNER="run_v1_two_barrier_dbtt_fatigue_map.py"
  else
    echo "error: could not find corrected two-barrier map runner." >&2
    echo "Set MAP_RUNNER=/path/to/run_v1_two_barrier_dbtt_fatigue_map_fixed.py" >&2
    exit 2
  fi
fi

if [[ ! -f "$CASE_TABLE" ]]; then
  echo "error: CASE_TABLE=$CASE_TABLE not found" >&2
  exit 2
fi

COMMON_PHYSICS=(
  --map-runner "$MAP_RUNNER"
  --case-table "$CASE_TABLE"
  --cycles-max "$CYCLES_MAX"
  --max-blocks "$MAX_BLOCKS"
  --criteria "$SENSITIVITY_RATE" "$PRIMARY_RATE"
  --DeltaK-min 0.025
  --DeltaK-max 20.0
  --DeltaK-seeds 0.05 0.10 0.20 0.40 0.80 1.60 3.20 6.40 12.80
  --threshold-abs-tol 0.05
  --threshold-rel-tol 0.03
  --max-refine-iters 10
  --resume
)

if [[ "$RUN_EXISTING_ANALYSIS" == "1" ]]; then
  if [[ -f "$EXISTING_MAP/fatigue_paris_points.csv" && -f "$EXISTING_MAP/monotonic_DBTT_points.csv" ]]; then
    echo "=== Re-analyzing existing corrected map with rate-defined thresholds ==="
    "$PYTHON" "$HERE/analyze_rate_defined_thresholds.py" \
      --map-dir "$EXISTING_MAP" \
      --out "$OUTROOT/existing_map_rate_analysis" \
      --criteria "$SENSITIVITY_RATE" "$PRIMARY_RATE" \
      --primary-criterion "$PRIMARY_RATE"
  else
    echo "warning: existing map CSVs not found at $EXISTING_MAP; skipping existing-map analysis" >&2
  fi
fi

if [[ "$RUN_CORE_ADAPTIVE" == "1" ]]; then
  echo "=== Adaptive core Kc(T) / DeltaK_th(T) study ==="
  # Six canonical cases; common Se=-40 and cleavage entropy {-5,0,+5}.
  "$PYTHON" "$HERE/run_adaptive_two_barrier_threshold_study.py" \
    "${COMMON_PHYSICS[@]}" \
    --out "$OUTROOT/core_adaptive" \
    --temperatures $CORE_TEMPS \
    --emission-entropies-kB -40 \
    --cleavage-entropies-kB -5 0 5

  "$PYTHON" "$HERE/analyze_rate_defined_thresholds.py" \
    --fatigue-csv "$OUTROOT/core_adaptive/fatigue_adaptive_points.csv" \
    --monotonic-csv "$OUTROOT/core_adaptive/monotonic_adaptive_points.csv" \
    --out "$OUTROOT/core_adaptive_analysis" \
    --criteria "$SENSITIVITY_RATE" "$PRIMARY_RATE" \
    --primary-criterion "$PRIMARY_RATE"
fi

if [[ "$RUN_EMISSION_STUDY" == "1" ]]; then
  echo "=== Emission-entropy fatigue study: cleavage entropy fixed at 0 kB ==="
  "$PYTHON" "$HERE/run_adaptive_two_barrier_threshold_study.py" \
    "${COMMON_PHYSICS[@]}" \
    --out "$OUTROOT/emission_adaptive" \
    --temperatures $EMISSION_TEMPS \
    --emission-entropies-kB $EMISSION_ENTROPIES \
    --cleavage-entropies-kB 0

  "$PYTHON" "$HERE/analyze_rate_defined_thresholds.py" \
    --fatigue-csv "$OUTROOT/emission_adaptive/fatigue_adaptive_points.csv" \
    --monotonic-csv "$OUTROOT/emission_adaptive/monotonic_adaptive_points.csv" \
    --out "$OUTROOT/emission_adaptive_analysis" \
    --criteria "$SENSITIVITY_RATE" "$PRIMARY_RATE" \
    --primary-criterion "$PRIMARY_RATE"
fi

if [[ "$RUN_STRENGTH" == "1" ]]; then
  echo "=== Arrhenius temperature-strength-anomaly study ==="
  "$PYTHON" "$HERE/run_arrhenius_strength_anomaly.py" \
    --map-runner "$MAP_RUNNER" \
    --case-table "$CASE_TABLE" \
    --out "$OUTROOT/strength_anomaly" \
    --emission-entropies-kB $EMISSION_ENTROPIES \
    --temperatures $STRENGTH_TEMPS \
    --strain-rates $STRENGTH_RATES \
    --strain-per-event "${STRAIN_PER_EVENT:-1.0}" \
    --sigma-max-GPa "${STRENGTH_SIGMA_MAX_GPA:-250}"
fi

if [[ "$RUN_CORRELATION" == "1" ]]; then
  LINK="$OUTROOT/emission_adaptive_analysis/rate_defined_DBTT_fatigue_link.csv"
  STRENGTH="$OUTROOT/strength_anomaly/strength_anomaly_metrics.csv"
  if [[ -f "$LINK" && -f "$STRENGTH" ]]; then
    echo "=== Strength-anomaly / excess-fatigue correlation analysis ==="
    "$PYTHON" "$HERE/analyze_strength_fatigue_correlation.py" \
      --link-csv "$LINK" \
      --strength-metrics "$STRENGTH" \
      --out "$OUTROOT/strength_fatigue_correlation" \
      --criterion "$PRIMARY_RATE"
  else
    echo "warning: correlation inputs missing; expected $LINK and $STRENGTH" >&2
  fi
fi

echo
echo "Study complete. Output root: $OUTROOT"
echo "Primary figures:"
echo "  $OUTROOT/core_adaptive_analysis/Kc_vs_rate_threshold_trajectories.png"
echo "  $OUTROOT/core_adaptive_analysis/fatigue_excess_residuals.png"
echo "  $OUTROOT/strength_anomaly/strength_temperature_anomaly_curves.png"
echo "  $OUTROOT/strength_fatigue_correlation/strength_anomaly_vs_fatigue_excess.png"
