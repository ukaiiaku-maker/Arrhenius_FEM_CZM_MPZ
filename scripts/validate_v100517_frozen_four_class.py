#!/usr/bin/env python3
"""Fail-closed preflight for the standalone frozen four-class FEM registry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from arrhenius_fracture.frozen_four_class_registry_v100517 import (
    EXPECTED_OPTIONS,
    load_frozen_parameter_option,
    validate_frozen_registry,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    payload = validate_frozen_registry()
    records = []
    for option in EXPECTED_OPTIONS:
        _, audit = load_frozen_parameter_option(option)
        records.append(audit)
    payload["records"] = records

    out = args.out.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    print(
        json.dumps(
            {
                "n_options": payload["n_options"],
                "classes": payload["classes"],
                "registry_sha256": payload["registry_sha256"],
                "active_parameter_fingerprint_sha256": payload[
                    "active_parameter_fingerprint_sha256"
                ],
                "phase_field_solver_dependency": False,
                "out": str(out),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
