from pathlib import Path
import re


def _executable_shell_text(text: str) -> str:
    """Return shell content with full-line and trailing comments removed."""
    executable = []
    for raw in text.splitlines():
        stripped = raw.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        executable.append(raw.split("#", 1)[0])
    return "\n".join(executable)


def _has_option(text: str, option: str) -> bool:
    return re.search(rf"(?<!\S){re.escape(option)}(?:=|\s|$)", text) is not None


def test_v1002_smoke_runner_uses_supported_straight_mode_i_cli():
    root = Path(__file__).resolve().parents[1]
    full_text = (
        root / "run_v10_0_2_progressive_one_segment_smoke.sh"
    ).read_text()
    commands = _executable_shell_text(full_text)

    assert "--da-phys 5e-6" in commands
    assert not _has_option(commands, "--da-phys-um")
    assert not _has_option(commands, "--no-crystal-branch")
    assert not _has_option(commands, "--crystal-branch")
    assert _has_option(commands, "--crystal-aniso")
    assert _has_option(commands, "--crystal-compete")
    assert "--max-fronts 1" in commands
    assert "mode_i_first_passage_v10_0_2_progressive" in commands

    # The anisotropic-elastic/straight-path semantics are tested against the
    # transformed function and runtime payload in the dedicated adapter tests.
