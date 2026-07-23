#!/usr/bin/env python3
"""Materialize the exact 4096-row registry from its versioned gzip asset."""

from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path
import tempfile


EXPECTED_ARCHIVE_SHA256 = (
    "59681e6c03b5bafdc07adadc140c6b3216ead5f7e1f60122dc0b4302658e0816"
)
EXPECTED_REGISTRY_SHA256 = (
    "1633851df78f4848a897d87e5ee9679f8e708095c1895369609bac2ab7c78efe"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--archive",
        type=Path,
        default=Path(
            "candidates/v9_12_targeted_local_4096_registry.csv.gz"
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("candidates/v9_12_targeted_local_4096_registry.csv"),
    )
    return parser.parse_args()


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def materialize(archive: Path, output: Path) -> None:
    archive_sha = _sha256_path(archive)
    if archive_sha != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError(
            f"candidate archive hash mismatch: {archive_sha} != "
            f"{EXPECTED_ARCHIVE_SHA256}"
        )
    if output.exists():
        output_sha = _sha256_path(output)
        if output_sha != EXPECTED_REGISTRY_SHA256:
            raise RuntimeError(
                f"existing candidate registry hash mismatch: {output_sha} != "
                f"{EXPECTED_REGISTRY_SHA256}; refusing to overwrite it"
            )
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=output.parent,
            prefix=f".{output.name}.",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            with gzip.open(archive, "rb") as source:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    temporary.write(block)
        output_sha = _sha256_path(temporary_path)
        if output_sha != EXPECTED_REGISTRY_SHA256:
            raise RuntimeError(
                f"materialized registry hash mismatch: {output_sha} != "
                f"{EXPECTED_REGISTRY_SHA256}"
            )
        temporary_path.replace(output)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def main() -> int:
    args = parse_args()
    materialize(args.archive, args.out)
    print(
        "V913_CANDIDATE_REGISTRY_OK "
        f"rows=4096 sha256={EXPECTED_REGISTRY_SHA256} out={args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
