#!/usr/bin/env bash
set -euo pipefail

OUT=${OUT:-runs/mpz_v9_10_4_narrow_dbtt_first_passage_fixed_zero_v1}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100"}
CLEAVAGE_SLOPE_MODE=${CLEAVAGE_SLOPE_MODE:-fixed_zero}
RESTARTS=${RESTARTS:-4}
DE_MAXITER=${DE_MAXITER:-80}
DE_POPSIZE=${DE_POPSIZE:-8}
LOCAL_MAXITER=${LOCAL_MAXITER:-250}
SEED=${SEED:-9104001}

printf '%s\n' \
  "========================================================================" \
  "v9.10.4.5 crash-safe narrow-DBTT first-passage search" \
  "temperatures=$TEMPS" \
  "restarts=$RESTARTS de_maxiter=$DE_MAXITER de_popsize=$DE_POPSIZE" \
  "local_maxiter=$LOCAL_MAXITER cleavage_slope_mode=$CLEAVAGE_SLOPE_MODE" \
  "out=$OUT" \
  "========================================================================"

PYTHONUNBUFFERED=1 python optimize_mpz_v9_10_4_5_narrow_dbtt.py \
  --temperatures "$TEMPS" \
  --cleavage-slope-mode "$CLEAVAGE_SLOPE_MODE" \
  --restarts "$RESTARTS" \
  --de-maxiter "$DE_MAXITER" \
  --de-popsize "$DE_POPSIZE" \
  --local-maxiter "$LOCAL_MAXITER" \
  --seed "$SEED" \
  --target-extension-um 5 \
  --out "$OUT"
