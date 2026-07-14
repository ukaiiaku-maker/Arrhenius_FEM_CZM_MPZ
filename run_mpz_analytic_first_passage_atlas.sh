#!/usr/bin/env bash
set -euo pipefail

# v9.2 analytical virgin-tip atlas. This stage evaluates the exact EXP-floor
# cleavage and emission hazards without advancing the moving process zone.

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_BIN="${PYTHON_BIN:-python}"
SCRIPT="${SCRIPT:-build_mpz_analytic_first_passage_atlas.py}"
INITIAL="${INITIAL:-mpz_three_class_initial_guesses.csv}"
OUTROOT="${OUTROOT:-runs/mpz_v9_2_analytic_first_passage_atlas}"
SHAPE_FAMILIES="${SHAPE_FAMILIES:-ceramic weakT DBTT}"
TEMPERATURES="${TEMPERATURES:-300 400 500 600 700 800 900 1000 1100 1200}"
KDOT_VALUES="${KDOT_VALUES:-0.005 0.02}"
SAMPLES_PER_FAMILY="${SAMPLES_PER_FAMILY:-32768}"
DK="${DK:-0.10}"
REFINE_DK="${REFINE_DK:-0.01}"
KMAX="${KMAX:-80}"
TOP_PER_REGION="${TOP_PER_REGION:-20}"
SEED="${SEED:-92031}"
PROGRESS_EVERY="${PROGRESS_EVERY:-200}"
FIRST_PARAMETERS="${FIRST_PARAMETERS:-runs/mpz_v9_1_three_class_tuning/01_first_passage/mpz_three_class_parameters.csv}"

export PYTHONUNBUFFERED=1

stamp() { date '+%Y-%m-%d %H:%M:%S'; }

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: analytical atlas script not found: $SCRIPT" >&2
  exit 2
fi
if [[ ! -f "$INITIAL" ]]; then
  echo "ERROR: initial shape table not found: $INITIAL" >&2
  exit 2
fi

ANCHOR_TABLES="$INITIAL"
if [[ -f "$FIRST_PARAMETERS" ]]; then
  ANCHOR_TABLES="$ANCHOR_TABLES $FIRST_PARAMETERS"
fi

mkdir -p "$OUTROOT"

echo "[$(stamp)] MPZ v9.2 analytical first-passage atlas"
echo "[$(stamp)] environment=$CONDA_ENV"
echo "[$(stamp)] output=$OUTROOT"
echo "[$(stamp)] shape_families=$SHAPE_FAMILIES"
echo "[$(stamp)] temperatures=$TEMPERATURES"
echo "[$(stamp)] Kdot_values=$KDOT_VALUES"
echo "[$(stamp)] samples_per_family=$SAMPLES_PER_FAMILY dK=$DK refine_dK=$REFINE_DK Kmax=$KMAX"
echo "[$(stamp)] anchor_tables=$ANCHOR_TABLES"

ARGS=(
  "$SCRIPT"
  --initial "$INITIAL"
  --shape-families "$SHAPE_FAMILIES"
  --temperatures "$TEMPERATURES"
  --Kdot-values "$KDOT_VALUES"
  --samples-per-family "$SAMPLES_PER_FAMILY"
  --seed "$SEED"
  --dK "$DK"
  --refine-dK "$REFINE_DK"
  --Kmax "$KMAX"
  --top-per-region "$TOP_PER_REGION"
  --anchor-tables "$ANCHOR_TABLES"
  --progress-every "$PROGRESS_EVERY"
  --out "$OUTROOT"
)

if command -v conda >/dev/null 2>&1; then
  conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u "${ARGS[@]}"
else
  "$PYTHON_BIN" -u "${ARGS[@]}"
fi

echo "[$(stamp)] analytical atlas complete: $OUTROOT"
