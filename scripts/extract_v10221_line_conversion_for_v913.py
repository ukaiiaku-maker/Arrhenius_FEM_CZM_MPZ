#!/usr/bin/env python3
"""Copy the signed 2-D kernel activation-to-line normalization into v9.13."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def arrays_for_key(value: Any, key: str) -> list[np.ndarray]:
    found: list[np.ndarray] = []
    if isinstance(value, dict):
        for name, child in value.items():
            if name == key:
                array = np.asarray(child, dtype=float).reshape(-1)
                if array.size:
                    found.append(array)
            found.extend(arrays_for_key(child, key))
    elif isinstance(value, list):
        for child in value:
            found.extend(arrays_for_key(child, key))
    return found


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signed-kernel-family", required=True)
    parser.add_argument("--base-physics", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-systems", type=int, default=2)
    args = parser.parse_args()

    family_path = Path(args.signed_kernel_family)
    physics_path = Path(args.base_physics)
    family = json.loads(family_path.read_text())
    arrays = arrays_for_key(family, "activation_to_line_content")
    valid = [
        array
        for array in arrays
        if array.shape == (args.n_systems,)
        and np.all(np.isfinite(array))
        and np.all(array > 0.0)
    ]
    if not valid:
        raise RuntimeError(
            "signed kernel family has no positive activation_to_line_content "
            f"array of length {args.n_systems}: {family_path}"
        )

    reference = valid[0]
    inconsistent = [
        array
        for array in valid[1:]
        if not np.allclose(array, reference, rtol=1e-10, atol=0.0)
    ]
    if inconsistent:
        unique = sorted({tuple(float(x) for x in array) for array in valid})
        raise RuntimeError(
            "activation_to_line_content is state dependent in this family; "
            "a single 1-D normalization is not exact. Unique arrays: "
            f"{unique[:8]}"
        )

    payload = json.loads(physics_path.read_text())
    common = payload.setdefault("common_physics", payload)
    common["activation_to_line_content_per_system"] = [
        float(value) for value in reference
    ]
    payload["conversion_provenance"] = str(family_path.resolve())
    payload["conversion_key"] = "activation_to_line_content"
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "V913_2D_LINE_CONVERSION_OK "
        f"values={','.join(f'{x:.17g}' for x in reference)} out={out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
