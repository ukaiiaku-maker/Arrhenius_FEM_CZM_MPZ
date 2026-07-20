"""v10.0.5.13.5 long-corridor startup repair.

The 20 um smoke corridor admitted a quality-valid physical-refinement mesh, but
at the 100 um production target the inherited candidate-count window started at
three centers.  With a 330 um radial cloud around every center, the additional
overlapping clouds can create Delaunay slivers.  This point release broadens the
deterministic search to include the lower-count candidates instead of lowering
the production triangle-quality floor.

No FEM, CZM, barrier, MPZ, hazard, source, shielding, or crack-growth law is
changed.
"""
from __future__ import annotations

import math
import sys

from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import mode_i_first_passage_v10_0_5_13_4_barrier_only as _base

POINT_RELEASE = "10.0.5.13.5"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_tip_source_long_corridor_v10_0_5_13_5"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST


def expanded_candidate_counts_v1005135(length_um: float, max_gap_um: float) -> list[int]:
    """Search all deterministic center counts from two through the old upper bound."""
    base = max(
        2,
        int(math.ceil(max(float(length_um), 0.0) / max(float(max_gap_um), 1.0))) + 1,
    )
    return list(range(2, base + 4))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    saved = _v91853._candidate_counts
    _v91853._candidate_counts = expanded_candidate_counts_v1005135
    try:
        return _base.main(user_args)
    finally:
        _v91853._candidate_counts = saved


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "expanded_candidate_counts_v1005135",
    "main",
]
