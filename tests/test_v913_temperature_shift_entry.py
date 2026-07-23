from __future__ import annotations

from pathlib import Path
import subprocess
import sys


def test_temperature_shift_script_can_be_invoked_directly() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "scripts" / "run_v913_dbtt_temperature_shift.py"
    completed = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "--target-peak-temperatures-K" in completed.stdout
    assert "--identity-tolerance-MPa-sqrt-m" in completed.stdout
