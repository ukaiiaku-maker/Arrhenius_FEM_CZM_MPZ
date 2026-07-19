"""v9.12.1 Mode-I entry with authoritative emitted-line field annotations."""
from __future__ import annotations

from . import mode_i_first_passage_v9_11 as _base
from . import sharp_front as _sharp_front
from .field_snapshots_v9121 import render_field_snapshots_v9121


def main(argv=None):
    original = _sharp_front._render_field_snapshots
    _sharp_front._render_field_snapshots = render_field_snapshots_v9121
    try:
        return _base.main(argv)
    finally:
        _sharp_front._render_field_snapshots = original


if __name__ == "__main__":
    main()
