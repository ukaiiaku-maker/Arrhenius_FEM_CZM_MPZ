from __future__ import annotations

import pytest

from arrhenius_fracture.mode_i_first_passage_v10_0_5_9_production_j_probe import (
    _ensure_v911_probe_contract,
)


def _value(args: list[str], option: str, default: str | None = None) -> str | None:
    try:
        index = args.index(option)
    except ValueError:
        return default
    return args[index + 1] if index + 1 < len(args) else default


def test_v911_probe_contract_adds_required_crystal_competition():
    args = _ensure_v911_probe_contract(
        ["--mode", "2d", "--crystal-aniso", "--max-fronts", "1"]
    )
    assert "--crystal-aniso" in args
    assert "--crystal-compete" in args
    assert "--crystal-branch" not in args
    assert _value(args, "--max-fronts") == "1"


def test_v911_probe_contract_adds_anisotropy_when_omitted():
    args = _ensure_v911_probe_contract(["--mode", "2d", "--max-fronts", "1"])
    assert "--crystal-aniso" in args
    assert "--crystal-compete" in args


def test_v911_probe_contract_rejects_incompatible_modes():
    with pytest.raises(SystemExit, match="requires --crystal-aniso"):
        _ensure_v911_probe_contract(
            ["--mode", "2d", "--no-crystal-aniso", "--max-fronts", "1"]
        )
    with pytest.raises(SystemExit, match="branching disabled"):
        _ensure_v911_probe_contract(
            [
                "--mode",
                "2d",
                "--crystal-aniso",
                "--crystal-branch",
                "--max-fronts",
                "1",
            ]
        )
    with pytest.raises(SystemExit, match="--max-fronts 1"):
        _ensure_v911_probe_contract(
            ["--mode", "2d", "--crystal-aniso", "--max-fronts", "2"]
        )
