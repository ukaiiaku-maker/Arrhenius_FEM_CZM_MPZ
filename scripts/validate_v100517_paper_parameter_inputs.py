#!/usr/bin/env python3
"""Resolve and audit all v10.0.5.17 paper parameter inputs before simulation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from arrhenius_fracture.audited_pf_parameter_bridge_v100517 import (
    load_audited_parameter_option,
)

ENTRY_BY_OPTION = {
    "v913_paper_peak01_0242980_persistent_sites": "v10.2.25",
    "v913_paper_peak02_0127508_persistent_sites": "v10.2.25",
    "v913_paper_peak03_0115460_persistent_sites": "v10.2.25",
    "v913_paper_dbtt01_0202500_persistent_sites": "v10.2.25",
    "v913_paper_dbtt02_0088403_persistent_sites": "v10.2.25",
    "v913_paper_control01_0086420_persistent_sites": "v10.2.25",
    "v913_paper_weakT01_0257068_persistent_sites": "v10.2.26",
    "v913_paper_ceramic01_0189364_persistent_sites": "v10.2.26",
}


def validate(parameter_source_root: Path, options: list[str], out: Path) -> dict:
    if len(options) != len(set(options)):
        raise ValueError("parameter option list contains duplicates")
    records = []
    candidate_ids = set()
    fingerprints = set()
    for option in options:
        try:
            entry = ENTRY_BY_OPTION[option]
        except KeyError as exc:
            raise KeyError(f"unknown paper parameter option {option!r}") from exc
        candidate, audit = load_audited_parameter_option(
            parameter_source_root,
            entry,
            option,
        )
        if candidate.candidate_id in candidate_ids:
            raise ValueError(f"duplicate candidate id {candidate.candidate_id}")
        if audit["selected_row_sha256"] in fingerprints:
            raise ValueError(f"duplicate active parameter fingerprint for {option}")
        candidate_ids.add(candidate.candidate_id)
        fingerprints.add(audit["selected_row_sha256"])
        audit["runtime_parameter_source_root"] = str(
            parameter_source_root.resolve()
        )
        audit["external_PF_runtime_dependency"] = False
        audit["parameter_source_vendored_into_FEM_CZM_checkout"] = True
        records.append(audit)
    payload = {
        "schema": "v10.0.5.17_audited_paper_parameter_preflight",
        "parameter_source_root": str(parameter_source_root.resolve()),
        "external_PF_runtime_dependency": False,
        "parameter_source_vendored_into_FEM_CZM_checkout": True,
        "n_options": len(records),
        "all_candidate_ids_unique": True,
        "all_selected_row_fingerprints_unique": True,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
        "options": records,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return payload


def _source_root(args: argparse.Namespace) -> Path:
    direct = args.parameter_source_root
    legacy = args.pf_repo_root
    if direct is not None and legacy is not None:
        if direct.expanduser().resolve() != legacy.expanduser().resolve():
            raise SystemExit(
                "--parameter-source-root and legacy --pf-repo-root disagree"
            )
    selected = direct if direct is not None else legacy
    if selected is None:
        raise SystemExit("--parameter-source-root is required")
    return selected.expanduser().resolve()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parameter-source-root", type=Path, default=None)
    parser.add_argument(
        "--pf-repo-root",
        type=Path,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("options", nargs="+")
    args = parser.parse_args()
    payload = validate(
        _source_root(args),
        args.options,
        args.out.expanduser().resolve(),
    )
    print(
        json.dumps(
            {
                "n_options": payload["n_options"],
                "parameter_source_root": payload["parameter_source_root"],
                "external_PF_runtime_dependency": False,
                "out": str(args.out),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
