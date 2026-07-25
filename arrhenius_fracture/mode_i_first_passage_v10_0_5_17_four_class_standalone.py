"""Run the frozen four-class rows on validated v10.0.5.17 FEM/CZM mechanics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v10_0_5_16_stochastic_pf_parity as _base
from .frozen_four_class_registry_v100517 import (
    BRIDGE_SCHEMA,
    load_frozen_parameter_option,
)
from .persistent_site_registry_v100514 import ROWS

POINT_RELEASE = "10.0.5.17-frozen-four-class"
MODEL_ID = "FEM_CZM_full_2D_frozen_four_class_stochastic_v10_0_5_17"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_17.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_17.json"
TRANSFER_MANIFEST = "frozen_four_class_parameter_transfer_v10_0_5_17.json"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "--parameter-entry",
        choices=("v10.0.5.17-frozen-four-class",),
        required=True,
    )
    parser.add_argument("--parameter-option", required=True)
    return parser


def _replace_option(argv: list[str], name: str, value: str) -> None:
    while name in argv:
        index = argv.index(name)
        del argv[index : min(index + 2, len(argv))]
    argv.extend([name, str(value)])


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _rewrite_output_metadata(out: Path | None, audit: dict[str, Any]) -> None:
    if out is None:
        return
    out.mkdir(parents=True, exist_ok=True)
    (out / TRANSFER_MANIFEST).write_text(
        json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n"
    )

    old_manifest = out / _base.PRODUCTION_MANIFEST
    if old_manifest.is_file():
        payload = json.loads(old_manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_17_frozen_four_class"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        payload["frozen_four_class_parameter_transfer"] = audit
        payload["parameter_entry"] = "v10.0.5.17-frozen-four-class"
        payload["parameter_option"] = audit["parameter_option"]
        payload["candidate_id"] = audit["candidate_id"]
        payload["material_class"] = audit.get("material_class")
        payload["paper_role"] = audit.get("paper_role")
        persistent = dict(payload.get("persistent_site_option", {}) or {})
        persistent.update(
            {
                "parameter_source": audit["registry_path"],
                "option_key": audit["parameter_option"],
                "candidate_id": audit["candidate_id"],
                "material_class": audit.get("material_class"),
                "paper_role": audit.get("paper_role"),
                "selected_row_sha256": audit["selected_row_sha256"],
                "active_parameter_fingerprint_sha256": audit[
                    "active_parameter_fingerprint_sha256"
                ],
                "selected_through_frozen_FEM_registry": True,
            }
        )
        payload["persistent_site_option"] = persistent
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(
            {
                "frozen_parameter_bridge_schema": BRIDGE_SCHEMA,
                "frozen_parameter_entry": "v10.0.5.17-frozen-four-class",
                "frozen_parameter_option": audit["parameter_option"],
                "frozen_candidate_id": audit["candidate_id"],
                "frozen_registry_sha256": audit["registry_sha256"],
                "frozen_active_parameter_fingerprint_sha256": audit[
                    "active_parameter_fingerprint_sha256"
                ],
                "phase_field_solver_dependency": False,
                "phase_field_result_comparison_enabled": False,
                "parameter_values_manually_reconstructed": False,
                "constitutive_parameters_changed_from_selected_registry_row": False,
                "mechanics_changed_for_parameter_installation": False,
                "source_closure_changed": False,
                "stochastic_cleavage_law_changed": False,
                "persistent_sites": True,
                "finite_source_inventory": False,
                "source_refresh": False,
                "explicit_recovery": False,
            }
        )
        payload["physics_contract"] = physics
        (out / PRODUCTION_MANIFEST).write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )
        old_manifest.unlink()

    old_selection = out / _base.SELECTION_MANIFEST
    if old_selection.is_file():
        selection = json.loads(old_selection.read_text())
        selection["schema"] = MODEL_ID
        selection["point_release"] = POINT_RELEASE
        selection["frozen_four_class_parameter_transfer"] = audit
        selection["parameter_entry"] = "v10.0.5.17-frozen-four-class"
        selection["parameter_option"] = audit["parameter_option"]
        selection["candidate_id"] = audit["candidate_id"]
        selection["material_class"] = audit.get("material_class")
        selection["paper_role"] = audit.get("paper_role")
        selection["phase_field_solver_dependency"] = False
        selection["phase_field_result_comparison_enabled"] = False
        (out / SELECTION_MANIFEST).write_text(
            json.dumps(selection, indent=2, sort_keys=True, default=str) + "\n"
        )
        old_selection.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    wrapper, remaining = _parser().parse_known_args(user_args)
    candidate, audit = load_frozen_parameter_option(wrapper.parameter_option)
    out = _out_path(remaining)
    _replace_option(remaining, "--persistent-site-option", candidate.option_key)

    existed = candidate.option_key in ROWS
    previous = ROWS.get(candidate.option_key)
    if existed and previous != candidate:
        raise RuntimeError(
            f"runtime registry already contains a different row for {candidate.option_key}"
        )
    ROWS[candidate.option_key] = candidate
    try:
        return _base.main(remaining)
    finally:
        _rewrite_output_metadata(out, audit)
        if existed:
            ROWS[candidate.option_key] = previous
        else:
            ROWS.pop(candidate.option_key, None)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
