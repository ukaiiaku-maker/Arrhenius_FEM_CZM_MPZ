"""v9.18.5.1 startup-safe committed-target horizon.

v9.18.5 introduced a dynamic step-horizon proxy so the accepted-step loop exits
immediately after the committed physical target.  The proxy implemented common
numeric operators but was not an actual ``int``.  Initialization paths in the
2-D driver and downstream configuration code are entitled to require strict
integer/index semantics, so the run could fail before the first FEM event.

This patch changes only that representation: the horizon is now an ``int``
subclass.  It behaves exactly like the requested nominal step count for range,
indexing, formatting, NumPy, JSON, and arithmetic, while its reflected
comparison ``step < horizon`` returns False once the controller requests the
post-commit stop.

No geometry, hazard, cohesive, MPZ, wake, shielding, or material law changes.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v9_18_5 as _v9185


class SafeDynamicStepHorizon(int):
    """Integer step horizon with one dynamic reflected comparison.

    Python gives the right-hand subclass reflected comparison priority for
    ``plain_int < SafeDynamicStepHorizon``.  All other integer operations are
    inherited directly from ``int`` and therefore retain strict index/numeric
    compatibility throughout solver initialization.
    """

    def __new__(cls, value: int, controller: Any):
        obj = int.__new__(cls, int(value))
        obj.controller = controller
        return obj

    def _running(self) -> bool:
        return not bool(getattr(self.controller, "v9185_stop_requested", False))

    def __gt__(self, other):
        # Reflected operation for ``accepted_step < horizon``.
        return self._running() and int(other) < int(self)

    def __ge__(self, other):
        return self._running() and int(other) <= int(self)



def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original = _v9185._DynamicStepHorizon
    _v9185._DynamicStepHorizon = SafeDynamicStepHorizon
    try:
        results = _v9185.main(user_args)
    finally:
        _v9185._DynamicStepHorizon = original

    out_value = _v9185._option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": "safe_target_stop_horizon_v91851_v1",
            "nominal_step_horizon_is_int_subclass": True,
            "strict_integer_index_semantics_preserved": True,
            "dynamic_post_commit_loop_comparison_enabled": True,
            "constitutive_physics_changed": False,
        }
        (out / "safe_target_stop_horizon_v91851.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
