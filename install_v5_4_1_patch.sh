#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-.}"
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"

if [[ ! -d "$ROOT/arrhenius_fracture" ]]; then
  echo "ERROR: $ROOT does not look like the Fatigue-PF project root (arrhenius_fracture/ missing)." >&2
  exit 2
fi

cp -f "$SCRIPT_DIR/run_sn_v1_representative_exp_floor_map.py" "$ROOT/"
cp -f "$SCRIPT_DIR/run_sn_v1_representative_exp_floor_map.sh" "$ROOT/"
cp -f "$SCRIPT_DIR/arrhenius_fracture/sn_arrhenius_chain.py" "$ROOT/arrhenius_fracture/"
cp -f "$SCRIPT_DIR/arrhenius_fracture/sn_v1_arrhenius_batch.py" "$ROOT/arrhenius_fracture/"
cp -f "$SCRIPT_DIR/arrhenius_fracture/sn_v1.py" "$ROOT/arrhenius_fracture/"
cp -f "$SCRIPT_DIR/arrhenius_fracture/sn_v1_arrhenius.py" "$ROOT/arrhenius_fracture/"

cd "$ROOT"
python - <<'PY'
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_v1_arrhenius import SNCase
from arrhenius_fracture.sn_v1_arrhenius_batch import run_stress_grid
print("V5.4.1 import preflight OK")
PY

echo "Installed V5.4.1 files into: $(pwd)"
