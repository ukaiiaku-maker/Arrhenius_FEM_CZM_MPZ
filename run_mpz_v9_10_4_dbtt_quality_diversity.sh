#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
SEARCH_ROOT="${SEARCH_ROOT:-runs/mpz_v9_10_3_dbtt_targeted_global_search_v1}"
OUTROOT="${OUTROOT:-runs/mpz_v9_10_4_dbtt_quality_diversity_v1}"
COUNT="${COUNT:-10}"
QUALITY_RESERVE_FRACTION="${QUALITY_RESERVE_FRACTION:-0.30}"
QUALITY_WEIGHT="${QUALITY_WEIGHT:-0.35}"
PARAMETER_WEIGHT="${PARAMETER_WEIGHT:-0.45}"
RESPONSE_WEIGHT="${RESPONSE_WEIGHT:-0.55}"
POOL_FACTOR="${POOL_FACTOR:-12}"
FORCE="${FORCE:-0}"

args=(
  --search-root "$SEARCH_ROOT"
  --target-class DBTT
  --out "$OUTROOT"
  --count "$COUNT"
  --quality-reserve-fraction "$QUALITY_RESERVE_FRACTION"
  --quality-weight "$QUALITY_WEIGHT"
  --parameter-weight "$PARAMETER_WEIGHT"
  --response-weight "$RESPONSE_WEIGHT"
  --pool-factor "$POOL_FACTOR"
)
if [[ "$FORCE" == "1" ]]; then
  args+=(--force)
fi

export PYTHONUNBUFFERED=1
conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u \
  select_mpz_v9_10_4_quality_diversity.py \
  "${args[@]}"
