#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PFROOT=${PFROOT:-/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_22_physical_front_width_top5_dbtt_screen}
PARAMETER_SET=${PARAMETER_SET:-paper4}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_17_${PARAMETER_SET}_anisotropic_300_1200K_200um_stochastic_pf_parity_v1}
INCLUDE_CONTROL=${INCLUDE_CONTROL:-0}

PAPER4_OPTIONS=(
  v913_paper_peak01_0242980_persistent_sites
  v913_paper_dbtt01_0202500_persistent_sites
  v913_paper_weakT01_0257068_persistent_sites
  v913_paper_ceramic01_0189364_persistent_sites
)
FULL7_OPTIONS=(
  v913_paper_peak01_0242980_persistent_sites
  v913_paper_peak02_0127508_persistent_sites
  v913_paper_peak03_0115460_persistent_sites
  v913_paper_dbtt01_0202500_persistent_sites
  v913_paper_dbtt02_0088403_persistent_sites
  v913_paper_weakT01_0257068_persistent_sites
  v913_paper_ceramic01_0189364_persistent_sites
)
case "$PARAMETER_SET" in
  paper4) OPTIONS=("${PAPER4_OPTIONS[@]}") ;;
  full7) OPTIONS=("${FULL7_OPTIONS[@]}") ;;
  *) echo "ERROR: PARAMETER_SET must be paper4 or full7" >&2; exit 1 ;;
esac
if [[ "$INCLUDE_CONTROL" == "1" ]]; then
  OPTIONS+=(v913_paper_control01_0086420_persistent_sites)
elif [[ "$INCLUDE_CONTROL" != "0" ]]; then
  echo "ERROR: INCLUDE_CONTROL must be 0 or 1" >&2
  exit 1
fi

mkdir -p "$CAMPAIGN_ROOT"
python "$ROOT/scripts/validate_v100517_paper_parameter_inputs.py" \
  --pf-repo-root "$PFROOT" \
  --out "$CAMPAIGN_ROOT/audited_parameter_catalog.json" \
  "${OPTIONS[@]}"

export ROOT PFROOT CAMPAIGN_ROOT INCLUDE_CONTROL PARAMETER_SET
exec bash "$ROOT/run_v10_0_5_17_paper7_300_1200K_200um_stochastic_campaign.sh"
