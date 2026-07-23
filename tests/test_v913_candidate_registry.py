from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from scripts.materialize_v913_candidate_registry import (
    EXPECTED_REGISTRY_SHA256,
    materialize,
)


def _sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_versioned_candidate_registry_materializes_exactly(tmp_path: Path):
    archive = Path(
        "candidates/v9_12_targeted_local_4096_registry.csv.gz"
    )
    output = tmp_path / "registry.csv"
    materialize(archive, output)
    assert _sha256_path(output) == EXPECTED_REGISTRY_SHA256
    with output.open() as stream:
        assert sum(1 for _ in stream) == 4097
    materialize(archive, output)


def test_materializer_refuses_to_replace_wrong_existing_registry(
    tmp_path: Path,
):
    archive = Path(
        "candidates/v9_12_targeted_local_4096_registry.csv.gz"
    )
    output = tmp_path / "registry.csv"
    output.write_text("not the accepted registry\n")
    with pytest.raises(RuntimeError, match="refusing to overwrite"):
        materialize(archive, output)
