#!/usr/bin/env python3
"""Apply the MPZ v9.1 core delta to an extracted MPZ v9.0 package."""
from __future__ import annotations

import argparse
import base64
import hashlib
import io
import shutil
from pathlib import Path
from zipfile import ZipFile

PAYLOAD_SHA256 = "8a49221c62faa9101e7fdc50968ddcb99aed3d74a24fd2f92885772a55126102"
PART_GLOB = "payload/core_delta.part*"
EXPECTED_FILES = 18


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base",
        required=True,
        type=Path,
        help="Extracted Arrhenius_FEM_CZM_MPZ_v9_0 directory",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="New v9.1 directory; defaults to a sibling versioned folder",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Apply directly to --base instead of making a copy",
    )
    args = parser.parse_args()

    repo_dir = Path(__file__).resolve().parent
    part_paths = sorted(repo_dir.glob(PART_GLOB))
    if len(part_paths) != 4:
        raise SystemExit(
            f"Expected 4 payload parts under {repo_dir / 'payload'}, found {len(part_paths)}"
        )

    encoded = "".join(p.read_text(encoding="ascii").strip() for p in part_paths)
    payload = base64.b64decode(encoded, validate=True)
    digest = hashlib.sha256(payload).hexdigest()
    if digest != PAYLOAD_SHA256:
        raise SystemExit(
            f"Payload checksum mismatch: expected {PAYLOAD_SHA256}, obtained {digest}"
        )

    base = args.base.expanduser().resolve()
    if not base.is_dir():
        raise SystemExit(f"Base directory not found: {base}")
    if args.in_place and args.out is not None:
        raise SystemExit("Use either --in-place or --out, not both")

    if args.in_place:
        target = base
    else:
        target = (
            args.out.expanduser().resolve()
            if args.out
            else base.with_name("Arrhenius_FEM_CZM_MPZ_v9_1_three_class_tuning")
        )
        if target.exists():
            raise SystemExit(f"Output already exists: {target}")
        print(f"Copying v9.0 to {target} ...")
        shutil.copytree(base, target)

    with ZipFile(io.BytesIO(payload), "r") as zf:
        bad = zf.testzip()
        if bad is not None:
            raise SystemExit(f"Payload is corrupt at: {bad}")
        members = [m for m in zf.infolist() if not m.is_dir()]
        if len(members) != EXPECTED_FILES:
            raise SystemExit(
                f"Expected {EXPECTED_FILES} payload files, found {len(members)}"
            )
        for member in members:
            rel = Path(member.filename)
            if rel.is_absolute() or ".." in rel.parts:
                raise SystemExit(f"Unsafe payload path: {member.filename}")
            dest = target / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(zf.read(member))

    (target / "APPLIED_MPZ_V9_1_DELTA.txt").write_text(
        "Applied MPZ v9.1 three-class tuning core delta.\n"
        f"Payload SHA-256: {PAYLOAD_SHA256}\n"
        f"Files applied: {len(members)}\n",
        encoding="utf-8",
    )

    print(f"Applied {len(members)} files to {target}")
    print(f"Payload SHA-256: {PAYLOAD_SHA256}")
    print("Next commands:")
    print(f"  cd {target}")
    print("  conda activate arrhenius-fem-czm")
    print("  python -m pip install -e .")
    print("  pytest -q tests/test_moving_process_zone.py tests/test_mpz_three_class_fit.py")
    print("  STAGE=smoke bash run_mpz_three_class_tuning.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
