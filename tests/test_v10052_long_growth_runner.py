from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "run_v10_0_5_2_DBTT_700K_100um_gate.sh"
RECOVERY = ROOT / "postprocess_v10_0_5_2_existing_100um_output.sh"


def _executable_text(path: Path) -> str:
    lines = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        lines.append(stripped)
    return "\n".join(lines)


def test_long_growth_runner_does_not_invoke_one_segment_audit():
    text = _executable_text(RUNNER)
    assert "audit_v10_0_3_progressive_integration.py" not in text
    assert "normalize_v10_0_3_1_reporting.py" in text
    assert "normalize_v10_0_5_1_slip_trace_reporting.py" in text
    assert "audit_v10_0_5_2_long_growth.py" in text
    assert text.index("normalize_v10_0_3_1_reporting.py") < text.index(
        "audit_v10_0_5_2_long_growth.py"
    )


def test_existing_output_recovery_is_postprocess_only():
    text = _executable_text(RECOVERY)
    assert "mode_i_first_passage" not in text
    assert "audit_v10_0_3_progressive_integration.py" not in text
    assert "run_completion_v10_0_5_2.json" in text
    assert "parallel_channel_diagnostics_v10_0_5_2.json" in text
    assert "normalize_v10_0_3_1_reporting.py" in text
    assert "normalize_v10_0_5_1_slip_trace_reporting.py" in text
    assert "audit_v10_0_5_2_long_growth.py" in text
