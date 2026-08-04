"""Fail-closed resolver for the frozen, immutable PF reference bundle.

The audited PF `steps_1000K.csv` and kernel `family.json` for the
Peak/1000K/theta=0/seed=8666 case became unavailable at their documented
path in the live PF reference repository during development (see
CLAUDE_PROGRESS.md's "EXTERNAL BLOCKER" and PF_REFERENCE_REGENERATION_CONTRACT.md).
This module defines where the *recovered, frozen* copies must live inside
this workspace, and refuses to proceed with anything else: no path here
ever points into a PF repository's own (mutable, externally-managed)
`runs/` tree.

Once the frozen bundle is in place, every production parity run must
resolve its PF inputs through this module rather than constructing paths
ad hoc, so a missing or silently-changed reference artifact fails the run
immediately instead of silently drifting onto different physics.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

MODEL_ID = "pf_theta0_frozen_reference_v10_0_5_18_4_0"

DEFAULT_REFERENCE_DIR = "reference_inputs/pf_v10_4_1_theta0_peak1000K_seed8666"

EXPECTED_STEPS_CSV_SHA256 = (
    "666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c"
)
EXPECTED_KERNEL_SHA256 = (
    "a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a"
)


@dataclass(frozen=True)
class FrozenPFReference:
    reference_dir: str
    steps_csv_path: str
    steps_csv_sha256: str
    kernel_path: str
    kernel_sha256: str


def _sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _verify_one(path: Path, *, expected_sha256: str, label: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(
            f"frozen PF reference {label} not found at {path}. "
            "This resolver never falls back to a live PF repository's "
            "runs/ tree -- see PF_REFERENCE_REGENERATION_CONTRACT.md for "
            "how to obtain and freeze this file."
        )
    observed = _sha256_of(path)
    if observed != expected_sha256:
        raise ValueError(
            f"frozen PF reference {label} at {path} has SHA-256 {observed}, "
            f"expected {expected_sha256}. Refusing to proceed with an "
            "artifact that does not match the audited contract -- do not "
            "relax this check; obtain the correct artifact instead."
        )
    return observed


def resolve_frozen_pf_reference(
    reference_dir: str | Path | None = None,
) -> FrozenPFReference:
    """Resolve and verify the frozen PF steps CSV and kernel family.json.

    `reference_dir` must be an explicit local directory (defaults to
    `DEFAULT_REFERENCE_DIR` under the current workspace). It is never
    inferred from, or allowed to silently fall back to, any path inside a
    PF repository's own `runs/` directory -- that directory is externally
    mutable and is exactly what made the original artifacts disappear.

    Raises FileNotFoundError or ValueError (fail closed) if either file is
    missing or its SHA-256 does not match the audited contract.
    """
    root = Path(reference_dir if reference_dir is not None else DEFAULT_REFERENCE_DIR)
    root = root.expanduser().resolve()

    steps_csv = root / "steps_1000K.csv"
    kernel = root / "family.json"

    steps_sha = _verify_one(
        steps_csv, expected_sha256=EXPECTED_STEPS_CSV_SHA256, label="steps_1000K.csv"
    )
    kernel_sha = _verify_one(
        kernel, expected_sha256=EXPECTED_KERNEL_SHA256, label="family.json"
    )

    return FrozenPFReference(
        reference_dir=str(root),
        steps_csv_path=str(steps_csv),
        steps_csv_sha256=steps_sha,
        kernel_path=str(kernel),
        kernel_sha256=kernel_sha,
    )


__all__ = [
    "MODEL_ID",
    "DEFAULT_REFERENCE_DIR",
    "EXPECTED_STEPS_CSV_SHA256",
    "EXPECTED_KERNEL_SHA256",
    "FrozenPFReference",
    "resolve_frozen_pf_reference",
]
