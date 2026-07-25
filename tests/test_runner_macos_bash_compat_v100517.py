from pathlib import Path
import subprocess

import numpy as np

from arrhenius_fracture.crystal import bcc_cleavage_traces


RUNNER = Path(
    "run_v10_0_5_17_paper7_300_1200K_200um_stochastic_campaign.sh"
)
WRAPPER4 = Path("run_v10_0_5_17_paper4_validated_campaign.sh")


def test_runner_uses_nonempty_command_array_under_nounset():
    text = RUNNER.read_text()
    assert "PLOT_ARGS" not in text
    assert "local -a CMD" in text
    assert "CMD=(" in text
    assert 'CMD+=(--no-plots)' in text
    assert '"${CMD[@]}"' in text


def test_four_class_default_and_direction_switch_gate():
    text = RUNNER.read_text()
    assert "PARAMETER_SET=${PARAMETER_SET:-paper4}" in text
    assert "PAPER4_OPTIONS=(" in text
    assert "v913_paper_peak01_0242980_persistent_sites" in text
    assert "v913_paper_dbtt01_0202500_persistent_sites" in text
    assert "v913_paper_weakT01_0257068_persistent_sites" in text
    assert "v913_paper_ceramic01_0189364_persistent_sites" in text
    assert "CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-30}" in text
    assert "--plane-gate-global" in text
    assert '--min-global-forward "$MIN_GLOBAL_FORWARD"' in text
    assert "--max-fronts 1" in text


def test_global_gate_allows_both_forward_BCC_100_traces_at_30deg():
    planes = bcc_cleavage_traces(theta_deg=30.0, include_110=False)
    forward = np.array([1.0, 0.0])
    directions = []
    for plane in planes:
        t = np.asarray(plane["t"], float)
        if float(t @ forward) < float((-t) @ forward):
            t = -t
        directions.append(t)
    assert all(float(t @ forward) >= 0.2 for t in directions)
    # The same alternate trace would be rejected by the old local-tangent gate.
    assert abs(float(directions[0] @ directions[1])) < 0.2


def test_runner_generates_partial_and_strict_analysis():
    text = RUNNER.read_text()
    assert "--allow-partial" in text
    assert "--lock-file" in text
    assert "RUN_STATUS=$?" in text
    assert "scripts/analyze_v100517_paper_parameter_campaign.py" in text


def test_paper4_validated_wrapper_selects_exactly_four_primary_options():
    text = WRAPPER4.read_text()
    assert "PARAMETER_SET=paper4" in text
    primary_block = text.split("OPTIONS=(", 1)[1].split("\n)", 1)[0]
    primary_options = [
        line.strip()
        for line in primary_block.splitlines()
        if line.strip().startswith("v913_paper_")
    ]
    assert primary_options == [
        "v913_paper_peak01_0242980_persistent_sites",
        "v913_paper_dbtt01_0202500_persistent_sites",
        "v913_paper_weakT01_0257068_persistent_sites",
        "v913_paper_ceramic01_0189364_persistent_sites",
    ]
    # The rehardening control is available only behind INCLUDE_CONTROL=1.
    assert "OPTIONS+=(v913_paper_control01_0086420_persistent_sites)" in text


def test_nonempty_array_expansion_is_safe_with_nounset(tmp_path):
    script = tmp_path / "probe.sh"
    script.write_text(
        """#!/usr/bin/env bash
set -eu
probe() {
  local mode="$1"
  local -a CMD
  CMD=(python -c 'print(1)')
  if [[ "$mode" == "0" ]]; then
    CMD+=(--no-plots)
  fi
  printf '%s:%s\n' "${#CMD[@]}" "${CMD[*]}"
}
probe 1
probe 0
"""
    )
    result = subprocess.run(
        ["bash", str(script)], check=True, text=True, capture_output=True
    )
    assert result.stdout.splitlines() == [
        "3:python -c print(1)",
        "4:python -c print(1) --no-plots",
    ]
