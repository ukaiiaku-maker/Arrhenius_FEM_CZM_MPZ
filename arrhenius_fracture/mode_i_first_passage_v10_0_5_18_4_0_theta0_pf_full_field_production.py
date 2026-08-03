"""Official v10.0.5.18.4.0 theta-zero full-field production entry.

The underlying parity solver executes first. Its physical state and files are
then left unchanged while inherited compatibility labels are normalized to the
audited theta-zero PF v10.4.1 full-field contract.

The direct prescribed-geometry PF kernel family may encode its explicitly
disabled wake channel with an empty source grid.  During this production entry
only, an audited schema adapter supplies the legacy loader with a one-point,
exactly-zero source wake representation.  The original PF artifact SHA, active
kernel, interpolation states, and runtime zero-wake physics remain unchanged.
"""
from __future__ import annotations

from pathlib import Path
import sys

from . import (
    mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_parity as _base,
)
from .active_only_kernel_family_compat_v10051840 import (
    AUDIT_FILE as KERNEL_COMPAT_AUDIT_FILE,
    installed_active_only_kernel_family_compat,
)
from .theta0_pf_full_field_metadata_v10051840 import (
    AUDIT_FILE as METADATA_AUDIT_FILE,
    normalize_theta0_full_field_outputs,
)

POINT_RELEASE = _base.POINT_RELEASE
MODEL_ID = (
    "FEM_CZM_theta0_PF_v10_4_1_full_field_production_v10_0_5_18_4_0"
)
AUDIT_FILE = _base.AUDIT_FILE
SEMANTIC_BULK_MODE = _base.SEMANTIC_BULK_MODE
SOLVER_BULK_MODE = _base.SOLVER_BULK_MODE


def _out_path(argv: list[str]) -> Path | None:
    raw = _base._control._option_value(argv, "--out")
    if raw is None:
        return None
    return Path(raw).expanduser().resolve()


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(args)
    if out is None:
        return _base.main(args)
    out.mkdir(parents=True, exist_ok=True)
    try:
        with installed_active_only_kernel_family_compat(out):
            return _base.main(args)
    finally:
        if out.exists():
            normalize_theta0_full_field_outputs(out)


if __name__ == "__main__":
    main()


__all__ = [
    "AUDIT_FILE",
    "KERNEL_COMPAT_AUDIT_FILE",
    "METADATA_AUDIT_FILE",
    "MODEL_ID",
    "POINT_RELEASE",
    "SEMANTIC_BULK_MODE",
    "SOLVER_BULK_MODE",
    "main",
]
