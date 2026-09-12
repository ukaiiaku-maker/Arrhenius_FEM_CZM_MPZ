#!/usr/bin/env python3
"""Export retained V5 evidence into the bounded one-dimensional V3 schema."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reduced_fracture_v3.oracle import export_v5_oracle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--static-evidence", type=Path, required=True)
    parser.add_argument("--checkpoint-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export_v5_oracle(
        static_evidence=args.static_evidence,
        checkpoint_directory=args.checkpoint_directory,
        output=args.output,
    )
    print(json.dumps(result["readiness"], sort_keys=True))


if __name__ == "__main__":
    main()
