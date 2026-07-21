#!/usr/bin/env python3
"""Convert a neutral/state-disabled 2-D trajectory into a v9.12 protocol CSV."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--extension-column", required=True)
    parser.add_argument("--K-column", required=True)
    parser.add_argument("--time-column")
    parser.add_argument("--default-duration-s", type=float, default=8.4)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with Path(args.input_csv).open(newline="") as fp:
        rows = list(csv.DictReader(fp))
    if len(rows) < 2:
        raise RuntimeError("input trajectory must contain at least two rows")
    points = []
    for row in rows:
        ext = float(row[args.extension_column])
        K = float(row[args.K_column])
        time = float(row[args.time_column]) if args.time_column else None
        points.append((ext, K, time))
    points.sort(key=lambda item: item[0])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fp:
        fields = [
            "extension_start_um",
            "extension_end_um",
            "K_start_MPa_sqrt_m",
            "K_end_MPa_sqrt_m",
            "duration_s",
        ]
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        for left, right in zip(points[:-1], points[1:]):
            if right[0] <= left[0]:
                continue
            if args.time_column:
                duration = right[2] - left[2]
                if duration <= 0.0:
                    raise RuntimeError("time column is not strictly increasing")
            else:
                duration = args.default_duration_s
            writer.writerow(
                {
                    "extension_start_um": f"{left[0]:.17g}",
                    "extension_end_um": f"{right[0]:.17g}",
                    "K_start_MPa_sqrt_m": f"{left[1]:.17g}",
                    "K_end_MPa_sqrt_m": f"{right[1]:.17g}",
                    "duration_s": f"{duration:.17g}",
                }
            )
    print(f"PROTOCOL_WRITTEN points={len(points)} out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
