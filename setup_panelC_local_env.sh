#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
ENV_PREFIX="${PANELC_ENV_PREFIX:-$PROJECT_ROOT/.conda-envs/panelC-sn}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda was not found on PATH." >&2
  exit 2
fi

if [[ -x "$ENV_PREFIX/bin/python" ]]; then
  echo "Panel C local environment already exists: $ENV_PREFIX"
else
  mkdir -p "$(dirname "$ENV_PREFIX")"
  conda create -y \
    --prefix "$ENV_PREFIX" \
    -c conda-forge \
    --strict-channel-priority \
    python=3.13 \
    'numpy<3' \
    'scipy<2' \
    'pandas<4' \
    'matplotlib<4'
fi

conda list --prefix "$ENV_PREFIX" --explicit \
  > "$PROJECT_ROOT/panelC_local_environment_explicit.txt"

echo "Local environment ready: $ENV_PREFIX"
"$ENV_PREFIX/bin/python" - <<'PY'
import sys, numpy, scipy, pandas, matplotlib
import scipy.sparse.linalg
print("python", sys.version.split()[0])
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("matplotlib", matplotlib.__version__)
print("scipy sparse import OK")
PY
