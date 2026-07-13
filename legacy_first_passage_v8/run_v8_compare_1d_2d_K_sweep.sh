#!/usr/bin/env bash
set -euo pipefail
OUT=${1:-runs/v8_compare_1d_2d_Ksweep_300K}
python run_v8_compare_1d_2d_K_sweep.py --out "$OUT" "$@"
