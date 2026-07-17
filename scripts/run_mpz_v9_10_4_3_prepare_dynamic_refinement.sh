#!/usr/bin/env bash
set -euo pipefail

INPUT_MANIFEST=${INPUT_MANIFEST:?Set INPUT_MANIFEST to the coarse first-passage shortlist}
OUT=${OUT:-runs/mpz_v9_10_4_3_dynamic_refinement_manifest_v1.csv}
COARSE_TEMPS=${COARSE_TEMPS:-"300 400 500 600 700 800 900 1000 1100"}
REFINEMENT_POINTS=${REFINEMENT_POINTS:-4}
SHELF_ANCHOR_COUNT=${SHELF_ANCHOR_COUNT:-2}
MAX_CANDIDATES=${MAX_CANDIDATES:-24}

PYTHONUNBUFFERED=1 python prepare_mpz_v9_10_4_3_dynamic_refinement.py \
  --input-manifest "$INPUT_MANIFEST" \
  --coarse-temperatures "$COARSE_TEMPS" \
  --refinement-points "$REFINEMENT_POINTS" \
  --shelf-anchor-count "$SHELF_ANCHOR_COUNT" \
  --max-candidates "$MAX_CANDIDATES" \
  --out "$OUT"
