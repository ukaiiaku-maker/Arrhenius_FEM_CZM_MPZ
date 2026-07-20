"""v10.0.5.13.3 tip/source-only barrier-transfer point release.

This wrapper changes no FEM, CZM, MPZ, barrier, source, shielding, transport, or
crack-advance implementation.  It records that the campaign must use the
existing ``tip_only`` v9.11 coupling: elastic continuum bulk plus the moving
crack-tip MPZ state.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_13_2_barrier_only as _base
from .mode_i_first_passage_v10_0 import _option_value

POINT_RELEASE = "10.0.5.13.3"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_tip_source_MPZ_v10_0_5_13_3"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST


def _update_manifest(out: Path, completed: bool) -> None:
    path = out / PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return
    payload.update(
        {
            "schema": "barrier_only_tip_source_manifest_v10_0_5_13_3",
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "plasticity_scope": {
                "bulk_plasticity_mode": "tip_only",
                "continuum_bulk_role": "elastic_fem_only",
                "moving_crack_tip_mpz_active": True,
                "uniform_bulk_mobile_retained_state_active": False,
                "candidate_source_or_shielding_closure_applied": False,
                "recorded_utc": datetime.now(timezone.utc).isoformat(),
            },
        }
    )
    if completed:
        payload["status"] = "complete"
        payload["run_completed_without_exception"] = True
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(user_args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.13.3 requires --out")
    mode = _option_value(user_args, "--bulk-plasticity-mode")
    if mode != "tip_only":
        raise SystemExit(
            "v10.0.5.13.3 requires --bulk-plasticity-mode tip_only; "
            f"received {mode!r}"
        )
    out = Path(out_value).resolve()
    try:
        result = _base.main(user_args)
        _update_manifest(out, completed=True)
        return result
    except BaseException:
        _update_manifest(out, completed=False)
        raise


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "MODEL_ID", "PRODUCTION_MANIFEST", "main"]
