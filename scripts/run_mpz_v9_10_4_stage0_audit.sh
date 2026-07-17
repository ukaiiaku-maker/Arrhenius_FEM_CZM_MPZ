#!/usr/bin/env bash
set -euo pipefail

MANIFEST=${MANIFEST:?Set MANIFEST to the historical DBTT spatial_promotion_manifest.csv}
OUT=${OUT:-runs/mpz_v9_10_4_3_current_dbtt_audit_v1}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100"}
TARGET_EXT_UM=${TARGET_EXT_UM:-5}
CANDIDATE_COUNT=${CANDIDATE_COUNT:-1}
HEARTBEAT_SECONDS=${HEARTBEAT_SECONDS:-30}

mkdir -p "$OUT"

echo "========================================================================"
echo "v9.10.4.3 DBTT historical-candidate audit"
echo "manifest=$MANIFEST"
echo "temperatures=$TEMPS"
echo "target_extension_um=$TARGET_EXT_UM candidate_count=$CANDIDATE_COUNT"
echo "heartbeat_seconds=$HEARTBEAT_SECONDS"
echo "out=$OUT"
echo "========================================================================"

PYTHONUNBUFFERED=1 python audit_mpz_v9_10_4_2_current_dbtt.py \
  --manifest "$MANIFEST" \
  --temperatures "$TEMPS" \
  --target-extension-um "$TARGET_EXT_UM" \
  --candidate-count "$CANDIDATE_COUNT" \
  --heartbeat-seconds "$HEARTBEAT_SECONDS" \
  --out "$OUT"
