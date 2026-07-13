#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV="${CONDA_ENV:-arrhenius-fem-czm}"
conda run -n "$ENV" python verify_mixed_mode_fem_czm_v8_1.py
