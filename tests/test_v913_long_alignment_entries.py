from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    "relative_path,required_option",
    [
        ("scripts/prepare_v913_long_peak25_campaign.py", "--source-ranking"),
        ("scripts/analyze_v913_long_peak_alignment.py", "--checkpoints-um"),
    ],
)
def test_entry_point_help(relative_path: str, required_option: str) -> None:
    root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [sys.executable, str(root / relative_path), "--help"],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert required_option in completed.stdout
