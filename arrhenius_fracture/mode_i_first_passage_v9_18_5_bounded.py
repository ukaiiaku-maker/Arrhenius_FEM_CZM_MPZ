"""Bounded wrapper for repeated v9.18.5 production-quality vetoes."""
from __future__ import annotations

import math
import os
import sys

import numpy as np

from . import mode_i_first_passage_v9_18_5 as _v9185


_actual_strict_advance = _v9185._strict_quality_advance


def _bounded_strict_quality_advance(self, *args, **kwargs):
    # v9.18.5.main assigns ``_original`` to whichever function is installed in
    # its module global. Forward that underlying backend entry point to the
    # actual strict-quality implementation before calling it.
    _actual_strict_advance._original = _bounded_strict_quality_advance._original
    result = _actual_strict_advance(self, *args, **kwargs)
    if bool(getattr(result, "inserted", False)):
        self._v9185_last_quality_veto_signature = None
        self._v9185_identical_quality_vetoes = 0
        return result

    reason = str(getattr(result, "reason", "unknown"))
    if not reason.startswith("v9185_quality_veto:"):
        return result

    front_id = int(kwargs.get("front_id", -1))
    p0 = np.asarray(kwargs.get("p0", [math.nan, math.nan]), float)
    p1 = np.asarray(kwargs.get("p1", [math.nan, math.nan]), float)
    scale = max(float(getattr(kwargs.get("mesh", None), "hbar_tip", 1.0e-12)), 1.0e-12)
    signature = (
        front_id,
        tuple(np.round(p0 / scale, 8)),
        tuple(np.round(p1 / scale, 8)),
        reason,
    )
    if signature == getattr(self, "_v9185_last_quality_veto_signature", None):
        count = int(getattr(self, "_v9185_identical_quality_vetoes", 0)) + 1
    else:
        count = 1
    self._v9185_last_quality_veto_signature = signature
    self._v9185_identical_quality_vetoes = count
    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        raise RuntimeError(
            "v9.18.5 repeated identical production-quality veto; physical renewal "
            f"remains unconsumed: front={front_id} count={count}/{limit} reason={reason}"
        )
    return result


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original = _v9185._strict_quality_advance
    _v9185._strict_quality_advance = _bounded_strict_quality_advance
    try:
        return _v9185.main(user_args)
    finally:
        _v9185._strict_quality_advance = original


if __name__ == "__main__":
    main()
