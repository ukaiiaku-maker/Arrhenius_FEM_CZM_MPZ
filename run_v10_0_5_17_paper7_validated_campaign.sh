#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PARAMETER_SOURCE_ROOT=${PARAMETER_SOURCE_ROOT:-$ROOT/runtime_inputs/v10_0_5_17_frozen_pf_inputs}
FAMILY_JSON=${FAMILY_JSON:-$PARAMETER_SOURCE_ROOT/signed_kernel/v10_2_14_active_only_campaign_family.json}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_17_paper7_300_1200K_200um_stochastic_pf_parity_v1}
INCLUDE_CONTROL=${INCLUDE_CONTROL:-0}

OPTIONS=(
  v913_paper_peak01_0242980_persistent_sites
  v913_paper_peak02_0127508_persistent_sites
  v913_paper_peak03_0115460_persistent_sites
  v913_paper_dbtt01_0202500_persistent_sites
  v913_paper_dbtt02_0088403_persistent_sites
  v913_paper_weakT01_0257068_persistent_sites
  v913_paper_ceramic01_0189364_persistent_sites
)
if [[ "$INCLUDE_CONTROL" == "1" ]]; then
  OPTIONS+=(v913_paper_control01_0086420_persistent_sites)
elif [[ "$INCLUDE_CONTROL" != "0" ]]; then
  echo "ERROR: INCLUDE_CONTROL must be 0 or 1" >&2
  exit 1
fi

mkdir -p "$CAMPAIGN_ROOT"
python "$ROOT/scripts/validate_v100517_paper_parameter_inputs.py" \
  --parameter-source-root "$PARAMETER_SOURCE_ROOT" \
  --out "$CAMPAIGN_ROOT/audited_parameter_catalog.json" \
  "${OPTIONS[@]}"

export ROOT PARAMETER_SOURCE_ROOT FAMILY_JSON CAMPAIGN_ROOT INCLUDE_CONTROL
exec bash "$ROOT/run_v10_0_5_17_paper7_300_1200K_200um_stochastic_campaign.sh"
