#!/usr/bin/env bash
set -euo pipefail

SOURCE=${SOURCE:?set SOURCE to the completed v9.10.4.8 dynamic_1d_all_candidates.csv}
OUT=${OUT:-runs/mpz_v9_10_4_9_six_candidate_ablation_v1}
TARGET_EXT_UM=${TARGET_EXT_UM:-5}
MODES=${MODES:-"full plasticity_off blunting_off backstress_off shielding_off"}
INCLUDE_BACKGROUND_FIELD_OFF=${INCLUDE_BACKGROUND_FIELD_OFF:-1}
EXPECTED_BRACKETS=${EXPECTED_BRACKETS:-6}

mkdir -p "$OUT"
MANIFEST="$OUT/six_candidate_ablation_manifest.csv"

echo "========================================================================"
echo "v9.10.4.9 bracket-balanced mechanism ablation"
echo "source=$SOURCE"
echo "out=$OUT"
echo "target_extension_um=$TARGET_EXT_UM"
echo "modes=$MODES"
echo "include_background_field_off=$INCLUDE_BACKGROUND_FIELD_OFF"
echo "expected_brackets=$EXPECTED_BRACKETS"
echo "========================================================================"

python select_mechanism_ablation_candidates_mpz_v9_10_4_9.py \
  --input "$SOURCE" \
  --out "$MANIFEST" \
  --per-bracket 1 \
  --expected-brackets "$EXPECTED_BRACKETS"

EXTRA_ARGS=()
if [[ "$INCLUDE_BACKGROUND_FIELD_OFF" == "1" ]]; then
  EXTRA_ARGS+=(--include-background-field-off)
fi

python evaluate_mechanism_ablation_mpz_v9_10_4_9.py \
  --manifest "$MANIFEST" \
  --out "$OUT" \
  --target-extension-um "$TARGET_EXT_UM" \
  --modes "$MODES" \
  --expected-candidates "$EXPECTED_BRACKETS" \
  "${EXTRA_ARGS[@]}"
