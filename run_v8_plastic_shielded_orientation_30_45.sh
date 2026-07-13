#!/usr/bin/env bash
set -euo pipefail

# Focused crystallographic fatigue-path study:
#   plastic_shielded_case64_M1
#   Kmax = 7.0 MPa sqrt(m)
#   theta = 30 and 45 deg
#
# The script preserves the validated v8 production driver and makes temporary
# patched copies that explicitly activate crystallographic path competition.
# Branching is enabled but not forced; it still requires the model's natural
# directional competition / hazard criteria.

DRIVER="${DRIVER:-run_v8_material_response_production_2d.sh}"

CASE_LABEL="${CASE_LABEL:-plastic_shielded_case64_M1}"
KMAX="${KMAX:-7.0}"

TARGET_EXT_UM="${TARGET_EXT_UM:-1470}"
BLOCKS="${BLOCKS:-25000}"
CYCLES_MAX="${CYCLES_MAX:-2e14}"
PRODUCTION_LEVEL="${PRODUCTION_LEVEL:-full}"
SNAPSHOTS="${SNAPSHOTS:-24}"
MAKE_2D_PLOTS="${MAKE_2D_PLOTS:-1}"

# Moderate crystallographic path anisotropy. Override from the command line
# if a sharper or weaker competition study is desired.
CLEAVE_GAMMA_ANISO="${CLEAVE_GAMMA_ANISO:-0.5}"
CRYSTAL_C44="${CRYSTAL_C44:-320e9}"

OUTBASE="${OUTBASE:-runs/v8_orientation_plastic_shielded_K7}"

if [[ ! -f "$DRIVER" ]]; then
  echo "ERROR: production driver not found: $DRIVER" >&2
  echo "Run this script from the Fatigue-PF project root." >&2
  exit 2
fi

mkdir -p "$OUTBASE"

TMPDIR_RUN="$(mktemp -d "${TMPDIR:-/tmp}/v8_orient.XXXXXX")"
trap 'rm -rf "$TMPDIR_RUN"' EXIT

make_oriented_driver() {
  local theta="$1"
  local src="$2"
  local dst="$3"

  python3 - "$src" "$dst" "$theta" "$CLEAVE_GAMMA_ANISO" "$CRYSTAL_C44" <<'PY'
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
theta = sys.argv[3]
gamma = sys.argv[4]
c44 = sys.argv[5]

text = src.read_text()

insert = (
    " --crystal-aniso"
    " --crystal-compete"
    " --crystal-branch"
    f" --crystal-theta-deg {theta}"
    f" --crystal-C44 {c44}"
    f" --cleave-gamma-aniso {gamma}"
)

# The validated v8 shell driver routes its 2-D calculations through commands
# containing '--mode 2d'. Inject into those commands without editing the source.
needle = "--mode 2d"
count = text.count(needle)

if count == 0:
    raise SystemExit(
        "ERROR: could not find '--mode 2d' in the production driver. "
        "The driver interface has changed, so crystal flags were NOT injected."
    )

patched = text.replace(needle, needle + insert)
dst.write_text(patched)
dst.chmod(0o755)

print(
    f"patched {count} occurrence(s) of '{needle}' for theta={theta} deg; "
    f"cleave_gamma_aniso={gamma}; C44={c44}"
)
PY
}

run_theta() {
  local theta="$1"
  local patched_driver="$TMPDIR_RUN/run_v8_theta${theta}.sh"
  local outroot="$OUTBASE/theta${theta}"

  make_oriented_driver "$theta" "$DRIVER" "$patched_driver"
  mkdir -p "$outroot"

  echo
  echo "======================================================================"
  echo "Orientation fatigue run"
  echo "  case              : $CASE_LABEL"
  echo "  Kmax              : $KMAX MPa sqrt(m)"
  echo "  crystal theta     : $theta deg"
  echo "  cleave anisotropy : $CLEAVE_GAMMA_ANISO"
  echo "  target extension  : $TARGET_EXT_UM um"
  echo "  blocks max        : $BLOCKS"
  echo "  cycles max        : $CYCLES_MAX"
  echo "  output root       : $outroot"
  echo "======================================================================"

  CASE_FILTER="$CASE_LABEL" \
  KLIST_OVERRIDE="$KMAX" \
  OUTROOT="$outroot" \
  PRODUCTION_LEVEL="$PRODUCTION_LEVEL" \
  TARGET_EXT_UM="$TARGET_EXT_UM" \
  BLOCKS="$BLOCKS" \
  CYCLES_MAX="$CYCLES_MAX" \
  SNAPSHOTS="$SNAPSHOTS" \
  MAKE_2D_PLOTS="$MAKE_2D_PLOTS" \
  bash "$patched_driver"
}

run_theta 30
run_theta 45

echo
echo "Completed both requested orientation runs."
echo "Results:"
echo "  $OUTBASE/theta30"
echo "  $OUTBASE/theta45"
