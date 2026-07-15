"""v9.13 Mode-I entry point with full-field and crack-tip zoom output."""
from __future__ import annotations

from . import mode_i_first_passage_v9_11 as _base
from . import sharp_front as _sharp_front
from .field_snapshots_v913 import render_field_snapshots_v913


def main(argv=None):
    original = _sharp_front._render_field_snapshots
    _sharp_front._render_field_snapshots = render_field_snapshots_v913
    try:
        return _base.main(argv)
    finally:
        _sharp_front._render_field_snapshots = original


if __name__ == "__main__":
    main()
