#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-}
TARGET_UM=${TARGET_UM:-100}
DA_UM=${DA_UM:-5}
MPZ_N_BINS=${MPZ_N_BINS:-200}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
LEGACY_RUN_COMMIT=${LEGACY_RUN_COMMIT:-4148a99}
LEGACY_GUARDED_RUNNER_ATTESTATION=${LEGACY_GUARDED_RUNNER_ATTESTATION:-0}

if [[ -z "$OUTROOT" ]]; then
  echo "ERROR: set OUTROOT to the completed v10.0.5.2 run directory"
  exit 1
fi
if [[ ! -d "$OUTROOT" ]]; then
  echo "ERROR: completed output directory not found: $OUTROOT"
  exit 1
fi
if [[ "$LEGACY_GUARDED_RUNNER_ATTESTATION" != "1" ]]; then
  echo "ERROR: this legacy output did not serialize the live MPZ bin count."
  echo "Set LEGACY_GUARDED_RUNNER_ATTESTATION=1 only for the output launched with"
  echo "run_v10_0_5_2_DBTT_700K_100um_gate.sh at commit $LEGACY_RUN_COMMIT."
  exit 1
fi
if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pip install -e . --no-deps

"$PYTHON_BIN" - "$OUTROOT" <<'PY'
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
required = [
    root / "run_completion_v10_0_5_2.json",
    root / "parallel_channel_diagnostics_v10_0_5_2.json",
    root / "mode_i_v10_0_5_2_results.json",
    root / "steps_0700K.csv",
]
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise SystemExit("missing completed v10.0.5.2 outputs: " + ", ".join(missing))
completion = json.loads(required[0].read_text())
if completion.get("status") != "complete" or completion.get("run_completed_without_exception") is not True:
    raise SystemExit("completion manifest does not certify a completed physical run")
print("Completed physical v10.0.5.2 run verified; no FEM solve will be launched.")
PY

# Reconstruct the authoritative outer v9.11 MPZ parser/factory binding. The
# inner run_args.json value is a compatibility shadow and is not the live state.
"$PYTHON_BIN" audit_v10_0_5_2_mpz_binding.py \
  "$OUTROOT" \
  --material DBTT \
  --expected-length-um "$MPZ_LENGTH_UM" \
  --expected-mpz-bins "$MPZ_N_BINS" \
  --run-commit "$LEGACY_RUN_COMMIT" \
  --guarded-runner-attested

# Do not invoke audit_v10_0_3_progressive_integration.py here. It is a
# deliberately one-segment smoke audit and rejects any valid multicommit run.
"$PYTHON_BIN" normalize_v10_0_3_1_reporting.py "$OUTROOT"
"$PYTHON_BIN" normalize_v10_0_5_1_slip_trace_reporting.py "$OUTROOT"
"$PYTHON_BIN" audit_v10_0_5_2_long_growth.py \
  "$OUTROOT" \
  --target-um "$TARGET_UM" \
  --expected-mpz-bins "$MPZ_N_BINS" \
  --da-um "$DA_UM"

cat <<EOF
V10.0.5.2 EXISTING 100 UM OUTPUT POSTPROCESSING PASSED
out=$OUTROOT
active_mpz_n_bins=$MPZ_N_BINS
active_mpz_length_um=$MPZ_LENGTH_UM
legacy_run_commit=$LEGACY_RUN_COMMIT
No FEM solve was launched.
The legacy one-segment audit was intentionally not used.
The inner run_args.json MPZ-bin value was retained only as a documented shadow default.
EOF
