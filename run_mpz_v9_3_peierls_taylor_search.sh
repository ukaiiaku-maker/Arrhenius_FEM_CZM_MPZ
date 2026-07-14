#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
ATLAS_SHORTLIST="${ATLAS_SHORTLIST:-runs/mpz_v9_2_analytic_first_passage_atlas/analytic_first_passage_atlas_shortlist_refined.csv}"
OUTROOT="${OUTROOT:-runs/mpz_v9_3_peierls_taylor_search}"
TRANSPORT_SAMPLES="${TRANSPORT_SAMPLES:-1024}"
INTRINSIC_TOP_PER_REGION="${INTRINSIC_TOP_PER_REGION:-5}"
TOP_PER_INTRINSIC="${TOP_PER_INTRINSIC:-3}"
TEMPERATURES="${TEMPERATURES:-300 700 900 1200}"
STRAIN_RATES="${STRAIN_RATES:-1e-5 1e-3}"
RHO_MIN="${RHO_MIN:-5e12}"
RHO_MAX="${RHO_MAX:-1e18}"
RHO_POINTS="${RHO_POINTS:-65}"
SEED="${SEED:-93017}"

export PYTHONUNBUFFERED=1
mkdir -p "$OUTROOT"

stamp() { date '+%Y-%m-%d %H:%M:%S'; }

echo "[$(stamp)] MPZ v9.3 emission-derived Peierls--Taylor search"
echo "[$(stamp)] atlas=$ATLAS_SHORTLIST"
echo "[$(stamp)] output=$OUTROOT"
echo "[$(stamp)] transport_samples=$TRANSPORT_SAMPLES intrinsic_top=$INTRINSIC_TOP_PER_REGION"
echo "[$(stamp)] temperatures=$TEMPERATURES strain_rates=$STRAIN_RATES"
echo "[$(stamp)] rho=$RHO_MIN..$RHO_MAX points=$RHO_POINTS"

ARGS=(
  search_mpz_peierls_taylor_parameters.py
  --atlas-shortlist "$ATLAS_SHORTLIST"
  --transport-samples "$TRANSPORT_SAMPLES"
  --intrinsic-top-per-region "$INTRINSIC_TOP_PER_REGION"
  --top-per-intrinsic "$TOP_PER_INTRINSIC"
  --temperatures "$TEMPERATURES"
  --strain-rates "$STRAIN_RATES"
  --rho-min "$RHO_MIN"
  --rho-max "$RHO_MAX"
  --rho-points "$RHO_POINTS"
  --seed "$SEED"
  --out "$OUTROOT"
)

if command -v conda >/dev/null 2>&1; then
  conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u "${ARGS[@]}"
else
  "$PYTHON_BIN" -u "${ARGS[@]}"
fi

echo "[$(stamp)] MPZ v9.3 Peierls--Taylor search complete"
