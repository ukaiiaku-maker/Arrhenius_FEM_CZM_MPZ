#!/usr/bin/env bash
set -euo pipefail

# Complete-ligament growth study: two conditions for each of the six canonical
# fatigue-barrier systems. Branching is not forced; the validated v8 multifront
# machinery remains responsible for natural branch birth and competition.

DRIVER="${DRIVER:-run_v8_material_response_production_2d.sh}"
OUTROOT="${OUTROOT:-runs/v8_material_response_complete_growth_2d}"

# Validated geometry: Lx=2.0 mm, initial notch a0=0.5 mm, completion near Lx-30 um.
TARGET_EXT_UM="${TARGET_EXT_UM:-1470}"

# Headroom above the ~1500 blocks used for ~150 um pilot runs.
BLOCKS="${BLOCKS:-25000}"
CYCLES_MAX="${CYCLES_MAX:-2e14}"
PRODUCTION_LEVEL="${PRODUCTION_LEVEL:-full}"
SNAPSHOTS="${SNAPSHOTS:-16}"
MAKE_2D_PLOTS="${MAKE_2D_PLOTS:-1}"

if [[ ! -f "$DRIVER" ]]; then
  echo "error: $DRIVER not found in $(pwd)" >&2
  echo "Run this wrapper from the project root containing the v8 production driver." >&2
  exit 2
fi

mkdir -p "$OUTROOT"

run_case() {
  local case_label="$1"
  local klist="$2"

  echo
  echo "======================================================================"
  echo "Complete-growth case: ${case_label}; Kmax=[${klist}]"
  echo "target extension=${TARGET_EXT_UM} um; blocks=${BLOCKS}; cycles_max=${CYCLES_MAX}"
  echo "======================================================================"

  CASE_FILTER="$case_label" \
  KLIST_OVERRIDE="$klist" \
  OUTROOT="$OUTROOT" \
  PRODUCTION_LEVEL="$PRODUCTION_LEVEL" \
  TARGET_EXT_UM="$TARGET_EXT_UM" \
  BLOCKS="$BLOCKS" \
  CYCLES_MAX="$CYCLES_MAX" \
  SNAPSHOTS="$SNAPSHOTS" \
  MAKE_2D_PLOTS="$MAKE_2D_PLOTS" \
  bash "$DRIVER"
}

# Moderate/stable-growth + faster/high-driving condition for each fatigue barrier.
run_case FCC_like_case29             "6.0 8.0"
run_case shifted_ductile_case64      "6.5 8.0"
run_case steep_cleavage_case35       "7.0 8.0"
run_case slow_threshold_case101      "7.5 8.0"
run_case higher_barrier_case171      "7.0 8.0"
run_case plastic_shielded_case64_M1  "7.0 8.0"

echo
echo "All 12 complete-growth runs have been passed to the v8 production driver."
echo "Output root: $OUTROOT"
