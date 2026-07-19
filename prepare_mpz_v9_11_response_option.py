#!/usr/bin/env python3
"""Materialize a runnable v9.11 parameter root for a named response option.

The option registry intentionally separates macroscopic response classes from
unique microscopic mechanisms.  In particular, a DBTT-like temperature trend
may be dominated by intrinsic first-passage kinetics, retained MPZ shielding,
or a mixture of both.  The generated root records the selected scientific role
in ``response_option_selection.json`` so archived runs remain interpretable.

This utility does not modify the canonical parameter root.  It copies the
canonical v9.11 root and replaces only the DBTT manifest when the selected
option is DBTT-based.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--option", required=True)
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("mpz_v9_11_response_options.json"),
    )
    parser.add_argument(
        "--dbtt-rows",
        type=Path,
        default=Path("mpz_v9_11_dbtt_option_rows.csv"),
    )
    parser.add_argument(
        "--canonical-root",
        type=Path,
        default=Path("mpz_v9_11_parameters"),
    )
    parser.add_argument("--outroot", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def load_registry(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise SystemExit(f"ERROR: option registry not found: {path}")
    payload = json.loads(path.read_text())
    if not isinstance(payload.get("options"), dict):
        raise SystemExit(f"ERROR: invalid option registry: {path}")
    return payload


def find_candidate_row(path: Path, candidate_id: str) -> tuple[list[str], dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: DBTT option row table not found: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: empty DBTT option row table: {path}")
        matches = [row for row in reader if row.get("candidate_id") == candidate_id]
    if len(matches) != 1:
        raise SystemExit(
            f"ERROR: expected one row for {candidate_id}; found {len(matches)} in {path}"
        )
    return list(reader.fieldnames), matches[0]


def require_canonical_root(root: Path) -> None:
    required = [
        root / "ceramic" / "spatial_promotion_manifest.csv",
        root / "weakT" / "spatial_promotion_manifest.csv",
        root / "DBTT" / "spatial_promotion_manifest.csv",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("ERROR: canonical parameter root is incomplete:\n  " + "\n  ".join(missing))


def main() -> int:
    args = parse_args()
    registry = load_registry(args.registry)
    options = registry["options"]
    if args.option not in options:
        available = ", ".join(sorted(options))
        raise SystemExit(f"ERROR: unknown option {args.option!r}; available: {available}")

    selection = dict(options[args.option])
    require_canonical_root(args.canonical_root)

    if args.outroot.exists():
        if not args.force:
            raise SystemExit(
                f"ERROR: destination already exists: {args.outroot}\n"
                "Use --force only when replacement is intentional."
            )
        shutil.rmtree(args.outroot)

    shutil.copytree(args.canonical_root, args.outroot)

    material_class = str(selection["material_class"])
    candidate_id = str(selection["candidate_id"])
    if material_class == "DBTT":
        fieldnames, row = find_candidate_row(args.dbtt_rows, candidate_id)
        row["target_class"] = "DBTT"
        if "selection_role" in row:
            row["selection_role"] = "primary"
        if "selection_reason" in row:
            row["selection_reason"] = f"v9_11_1_response_option:{args.option}"
        if "accepted_for_spatial_promotion" in row:
            row["accepted_for_spatial_promotion"] = "True"

        destination = args.outroot / "DBTT" / "spatial_promotion_manifest.csv"
        with destination.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({name: row.get(name, "") for name in fieldnames})
    elif material_class != "weakT":
        raise SystemExit(f"ERROR: unsupported material_class={material_class!r}")

    audit = {
        "schema_version": registry.get("schema_version", 1),
        "registry_version": registry.get("version"),
        "option": args.option,
        "selection": selection,
        "scientific_scope": registry.get("scientific_scope", {}),
        "canonical_root": str(args.canonical_root.resolve()),
        "dbtt_rows": str(args.dbtt_rows.resolve()),
        "outroot": str(args.outroot.resolve()),
    }
    (args.outroot / "response_option_selection.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    print("RESPONSE_OPTION_ROOT_READY")
    print(f"option={args.option}")
    print(f"candidate_id={candidate_id}")
    print(f"material_class={material_class}")
    print(f"role={selection['role']}")
    print(f"outroot={args.outroot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
