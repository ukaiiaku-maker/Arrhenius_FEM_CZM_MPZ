"""Audited v10.0.5.5 entry with complete cohesive cache signatures."""
from __future__ import annotations

from typing import Any

from . import mode_i_first_passage_v10_0_5_5_stochastic_vhcf as _base

POINT_RELEASE = _base.POINT_RELEASE
MODEL_ID = _base.MODEL_ID + "_audited"


_OLD_SIGNATURE = """        def _vhcf_cohesive_signature_v10055():
            if cohesive_network is None:
                return ()
            values = []
            for name in ('damage', 'd', 'state', 'normal_damage', 'tangential_damage'):
                value = getattr(cohesive_network, name, None)
                try:
                    arr = np.asarray(value, dtype=float)
                except Exception:
                    continue
                if arr.size > 0 and np.all(np.isfinite(arr)):
                    values.append((name, int(arr.size), float(np.sum(arr)), float(np.max(arr))))
            return tuple(values)
"""

_NEW_SIGNATURE = """        def _vhcf_cohesive_signature_v10055():
            if cohesive_network is None:
                return ()
            values = []
            for name in ('damage', 'd', 'state', 'normal_damage', 'tangential_damage'):
                value = getattr(cohesive_network, name, None)
                try:
                    arr = np.asarray(value, dtype=float)
                except Exception:
                    continue
                if arr.size > 0 and np.all(np.isfinite(arr)):
                    values.append((name, int(arr.size), float(np.sum(arr)), float(np.max(arr))))
            elements = list(getattr(cohesive_network, 'elements', []) or [])
            if elements:
                topology = []
                damage_sum = 0.0
                clock_sum = 0.0
                for elem in elements:
                    plus = tuple(int(x) for x in getattr(elem, 'plus_nodes', ()))
                    minus = tuple(int(x) for x in getattr(elem, 'minus_nodes', ()))
                    topology.append((plus, minus, int(getattr(elem, 'front_id', -1)),
                                     int(getattr(elem, 'event_index', -1))))
                    damage_sum += float(getattr(elem, 'damage', 0.0))
                    clock_sum += float(getattr(elem, 'clock', 0.0))
                values.append(('cohesive_elements', len(elements), tuple(topology),
                               float(damage_sum), float(clock_sum)))
            return tuple(values)
"""


def patch_run_2d_source_v10055_audited(source: str) -> str:
    patched = _ORIGINAL_PATCH(source)
    count = patched.count(_OLD_SIGNATURE)
    if count != 1:
        raise RuntimeError(
            "v10.0.5.5 audited cache signature expected one anchor; "
            f"found {count}"
        )
    return patched.replace(_OLD_SIGNATURE, _NEW_SIGNATURE)


_ORIGINAL_PATCH = _base.patch_run_2d_source_v10055


def validate_source_transform_v10055() -> dict[str, Any]:
    saved = _base.patch_run_2d_source_v10055
    _base.patch_run_2d_source_v10055 = patch_run_2d_source_v10055_audited
    try:
        result = dict(_base.validate_source_transform_v10055())
    finally:
        _base.patch_run_2d_source_v10055 = saved
    result.update(
        {
            "cohesive_element_cache_signature": True,
            "v10055_audited_entry": True,
        }
    )
    return result


def main(argv: list[str] | None = None):
    validate_source_transform_v10055()
    saved = _base.patch_run_2d_source_v10055
    _base.patch_run_2d_source_v10055 = patch_run_2d_source_v10055_audited
    try:
        return _base.main(argv)
    finally:
        _base.patch_run_2d_source_v10055 = saved


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "patch_run_2d_source_v10055_audited",
    "validate_source_transform_v10055",
    "main",
]
