#!/usr/bin/env bash
set -euo pipefail

OUT=${OUT:-runs/mpz_v9_10_4_7_analytical_downselect_v1}
SAMPLES=${SAMPLES:-8192}
SEED=${SEED:-9104701}
WORKERS=${WORKERS:-4}
BATCH_SIZE=${BATCH_SIZE:-128}
PER_BRACKET_KEEP=${PER_BRACKET_KEEP:-8}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100"}
CLEAVAGE_SLOPE_MODE=${CLEAVAGE_SLOPE_MODE:-fixed_zero}

mkdir -p "$OUT"

echo "========================================================================"
echo "v9.10.4.7 analytical DBTT down-selection"
echo "samples=$SAMPLES workers=$WORKERS batch_size=$BATCH_SIZE"
echo "per_bracket_keep=$PER_BRACKET_KEEP"
echo "temperatures=$TEMPS"
echo "cleavage_slope_mode=$CLEAVAGE_SLOPE_MODE"
echo "out=$OUT"
echo "========================================================================"

PYTHONUNBUFFERED=1 python analytical_downselect_mpz_v9_10_4_7.py \
  --samples "$SAMPLES" \
  --seed "$SEED" \
  --workers "$WORKERS" \
  --batch-size "$BATCH_SIZE" \
  --per-bracket-keep "$PER_BRACKET_KEEP" \
  --temperatures "$TEMPS" \
  --cleavage-slope-mode "$CLEAVAGE_SLOPE_MODE" \
  --out "$OUT"
