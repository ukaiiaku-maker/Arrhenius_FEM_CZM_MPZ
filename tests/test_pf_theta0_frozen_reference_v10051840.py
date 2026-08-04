from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from arrhenius_fracture import pf_theta0_frozen_reference_v10051840 as frozen_ref


def _write(path: Path, content: bytes) -> str:
    path.write_bytes(content)
    return hashlib.sha256(content).hexdigest()


def test_resolves_when_both_files_match_expected_hashes(tmp_path: Path, monkeypatch):
    csv_bytes = b"synthetic steps csv content"
    kernel_bytes = b"synthetic kernel content"
    monkeypatch.setattr(frozen_ref, "EXPECTED_STEPS_CSV_SHA256", hashlib.sha256(csv_bytes).hexdigest())
    monkeypatch.setattr(frozen_ref, "EXPECTED_KERNEL_SHA256", hashlib.sha256(kernel_bytes).hexdigest())
    (tmp_path / "steps_1000K.csv").write_bytes(csv_bytes)
    (tmp_path / "family.json").write_bytes(kernel_bytes)

    result = frozen_ref.resolve_frozen_pf_reference(tmp_path)
    assert result.reference_dir == str(tmp_path.resolve())
    assert result.steps_csv_sha256 == hashlib.sha256(csv_bytes).hexdigest()
    assert result.kernel_sha256 == hashlib.sha256(kernel_bytes).hexdigest()


def test_fails_closed_when_csv_missing(tmp_path: Path):
    (tmp_path / "family.json").write_bytes(b"x")
    with pytest.raises(FileNotFoundError, match="steps_1000K.csv"):
        frozen_ref.resolve_frozen_pf_reference(tmp_path)


def test_fails_closed_when_kernel_missing(tmp_path: Path, monkeypatch):
    csv_bytes = b"synthetic steps csv content"
    monkeypatch.setattr(frozen_ref, "EXPECTED_STEPS_CSV_SHA256", hashlib.sha256(csv_bytes).hexdigest())
    (tmp_path / "steps_1000K.csv").write_bytes(csv_bytes)
    with pytest.raises(FileNotFoundError, match="family.json"):
        frozen_ref.resolve_frozen_pf_reference(tmp_path)


def test_fails_closed_on_csv_hash_mismatch(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(frozen_ref, "EXPECTED_KERNEL_SHA256", hashlib.sha256(b"k").hexdigest())
    (tmp_path / "steps_1000K.csv").write_bytes(b"wrong content")
    (tmp_path / "family.json").write_bytes(b"k")
    with pytest.raises(ValueError, match="steps_1000K.csv"):
        frozen_ref.resolve_frozen_pf_reference(tmp_path)


def test_fails_closed_on_kernel_hash_mismatch(tmp_path: Path, monkeypatch):
    csv_bytes = b"synthetic steps csv content"
    monkeypatch.setattr(frozen_ref, "EXPECTED_STEPS_CSV_SHA256", hashlib.sha256(csv_bytes).hexdigest())
    (tmp_path / "steps_1000K.csv").write_bytes(csv_bytes)
    (tmp_path / "family.json").write_bytes(b"wrong kernel content")
    with pytest.raises(ValueError, match="family.json"):
        frozen_ref.resolve_frozen_pf_reference(tmp_path)


def test_never_falls_back_to_a_pf_repository_runs_path(tmp_path: Path):
    # An empty explicit directory must fail closed, never silently resolve
    # to some other conventional PF-repo-like path.
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        frozen_ref.resolve_frozen_pf_reference(empty_dir)


def test_default_reference_dir_is_local_to_this_workspace_not_a_pf_repo():
    assert "PF-fracture-fatigue" not in frozen_ref.DEFAULT_REFERENCE_DIR
    assert "runs" not in Path(frozen_ref.DEFAULT_REFERENCE_DIR).parts
    assert frozen_ref.DEFAULT_REFERENCE_DIR.startswith("reference_inputs/")


def test_expected_hashes_match_the_audited_contract():
    # Pin these literally so any accidental edit to the module is caught --
    # these two values are the entire point of this file and must never
    # drift silently.
    assert frozen_ref.EXPECTED_STEPS_CSV_SHA256 == (
        "666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c"
    )
    assert frozen_ref.EXPECTED_KERNEL_SHA256 == (
        "a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a"
    )
