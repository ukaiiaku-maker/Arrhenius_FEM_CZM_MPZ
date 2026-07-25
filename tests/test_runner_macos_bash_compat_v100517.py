from pathlib import Path
import subprocess


RUNNER = Path(
    "run_v10_0_5_17_paper7_300_1200K_200um_stochastic_campaign.sh"
)


def test_runner_uses_nonempty_command_array_under_nounset():
    text = RUNNER.read_text()
    assert "PLOT_ARGS" not in text
    assert "local -a CMD" in text
    assert "CMD=(" in text
    assert 'CMD+=(--no-plots)' in text
    assert '"${CMD[@]}"' in text


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
