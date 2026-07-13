#!/usr/bin/env python3
"""Patch sn_geometry.py for NumPy >=2.5, where np.cross no longer accepts 2-vectors.

Run from the Fatigue-PF repository root:
    python patch_sn_geometry_numpy25_cross2d.py

The script creates arrhenius_fracture/sn_geometry.py.bak before editing.
"""
from __future__ import annotations

from pathlib import Path
import shutil
import sys

path = Path("arrhenius_fracture/sn_geometry.py")
if not path.exists():
    sys.exit(f"ERROR: {path} not found. Run this script from the Fatigue-PF repository root.")

old = "    area2 = abs(np.cross(p2 - p1, p3 - p1))\n"
new = (
    "    v21 = p2 - p1\n"
    "    v31 = p3 - p1\n"
    "    # 2-D scalar cross product (twice the signed triangle area).\n"
    "    # np.cross on 2-vectors was deprecated in NumPy 2.0 and removed in 2.5.\n"
    "    area2 = abs(v21[0] * v31[1] - v21[1] * v31[0])\n"
)

text = path.read_text()
if new in text:
    print(f"Already patched: {path}")
    raise SystemExit(0)
if old not in text:
    sys.exit(
        "ERROR: expected np.cross line not found; file was not modified.\n"
        "Search manually with: grep -n 'np.cross' arrhenius_fracture/sn_geometry.py"
    )

backup = path.with_suffix(path.suffix + ".bak")
if not backup.exists():
    shutil.copy2(path, backup)
    print(f"backup: {backup}")

path.write_text(text.replace(old, new, 1))
print(f"patched: {path}")
print("next: python -m compileall -q arrhenius_fracture")
