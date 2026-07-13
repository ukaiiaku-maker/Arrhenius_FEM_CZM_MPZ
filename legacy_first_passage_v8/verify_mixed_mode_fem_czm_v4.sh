#!/usr/bin/env bash
set -euo pipefail
E="${CONDA_ENV:-arrhenius-fem-czm}"
conda run -n "$E" python -c 'import sys,numpy,scipy;print("python:",sys.executable);print("numpy:",numpy.__version__);print("scipy:",scipy.__version__)'
conda run -n "$E" python verify_mixed_mode_fem_czm_v4.py
