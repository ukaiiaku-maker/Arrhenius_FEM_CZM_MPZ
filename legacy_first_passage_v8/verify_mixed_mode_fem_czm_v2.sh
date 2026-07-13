#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)";cd "$ROOT"
CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
conda run -n "$CONDA_ENV" python verify_mixed_mode_fem_czm_v2.py
