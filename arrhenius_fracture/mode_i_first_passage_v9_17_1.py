"""v9.17.1 routing fix for the absolute-hazard/source-refresh controller.

v9.17 correctly requires one cleavage renewal per cohesive-opening trajectory,
but attempted to enforce that internal front configuration through an option
that is not exposed by the active v9.11 command-line parser.  This wrapper:

1. suppresses only the unsupported ``--max-advances-per-step`` injection; and
2. sets ``front.f.max_advances_per_step = 1`` when each MPZ engine is created.

No v9.17 constitutive law, hazard-clock rule, source-refresh rule, or audit
threshold is changed.
"""
from __future__ import annotations

import sys
from typing import Any, Callable

from . import mode_i_first_passage_v9_17 as _v917
from .mpz_front_engine_v911 import MovingProcessZone2DFrontEngine


def _one_fire_engine_initializer(original_init: Callable[..., Any]):
    """Return an engine initializer that enforces one renewal per solver step."""

    def wrapped(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        self.f.max_advances_per_step = 1.0
        self.v9171_one_fire_internal_routing = True

    return wrapped


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_inject = _v917._v916._inject_once
    original_init = MovingProcessZone2DFrontEngine.__init__

    def inject_without_unexposed_one_fire_option(
        args: list[str], name: str, value: str | None = None
    ) -> None:
        if name == "--max-advances-per-step":
            return
        original_inject(args, name, value)

    _v917._v916._inject_once = inject_without_unexposed_one_fire_option
    MovingProcessZone2DFrontEngine.__init__ = _one_fire_engine_initializer(original_init)
    try:
        return _v917.main(user_args)
    finally:
        _v917._v916._inject_once = original_inject
        MovingProcessZone2DFrontEngine.__init__ = original_init


if __name__ == "__main__":
    main()
