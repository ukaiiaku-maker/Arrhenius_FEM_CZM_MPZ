"""v9.14 direct Mode-I entry using the conservative adaptive-CZM event path.

The legacy direct Mode-I runner used the non-directional single-front branch in
``sharp_front``.  That branch inserted cohesive topology but did not execute the
mature ``_advance_polyline`` history-transfer path used by the anisotropic
multi-front solver.  v9.14 enables that path while constraining the admissible
cleavage inventory to one exactly forward Mode-I plane.  Branching remains off.

Consequences for every accepted cleavage event:

* the Arrhenius adaptive clock locates one event surface;
* one physical cohesive increment is inserted;
* local r/h adaptation may refine the tip patch;
* every new Gauss point inherits its parent element's history through
  ``elem_parent_map``;
* the explicit crack path and cohesive network remain authoritative.
"""
from __future__ import annotations

import sys
from typing import Any

import numpy as np

from . import crystal as _crystal
from . import mode_i_first_passage_v9_13 as _base


def _forward_mode_i_plane(theta_deg: float, *args: Any, **kwargs: Any):
    del theta_deg, args, kwargs
    return [{
        "name": "v914_forward_mode_I",
        "family": "mode_I_forward",
        "angle_deg": 0.0,
        "t": np.array([1.0, 0.0], dtype=float),
        "n": np.array([0.0, 1.0], dtype=float),
        "gamma_rel": 1.0,
    }]


def _inject_once(argv: list[str], name: str, value: str | None = None) -> None:
    if any(token == name or token.startswith(name + "=") for token in argv):
        return
    argv.append(name)
    if value is not None:
        argv.append(value)


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _inject_once(user_args, "--crystal-aniso")
    _inject_once(user_args, "--no-crystal-branch")
    _inject_once(user_args, "--crack-backend", "adaptive_czm")
    _inject_once(user_args, "--adaptive-events")
    _inject_once(user_args, "--crystal-compete", "0")
    _inject_once(user_args, "--min-global-forward", "0.999999")

    original = _crystal.bcc_cleavage_traces
    _crystal.bcc_cleavage_traces = _forward_mode_i_plane
    try:
        return _base.main(user_args)
    finally:
        _crystal.bcc_cleavage_traces = original


if __name__ == "__main__":
    main()
