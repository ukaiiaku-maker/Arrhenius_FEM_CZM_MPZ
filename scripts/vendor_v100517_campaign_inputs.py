#!/usr/bin/env python3
"""Vendor the exact frozen inputs used by the v10.0.5.17 paper7 campaign.

This utility does not modify any FEM, CZM, stochastic-event, or moving-process-zone
code. It reads the provenance written by the completed campaign, copies the
parameter registries and option maps from the exact recorded parameter-source
commit, and copies the independently recorded signed-kernel family byte-for-byte
from its recorded path after verifying its frozen SHA-256.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from arrhenius_fracture.audited_pf_parameter_bridge_v100517 import (
    load_audited_parameter_option,
)

BASELINE_COMMIT = "75fc7f1ab1a82f0bb631d859f791cda3709e60fd"
DEFAULT_DESTINATION = ROOT / "runtime_inputs" / "v10_0_5_17_frozen_pf_inputs"
CONFIGURATION_NAME = "campaign_configuration.txt"
MANIFEST_NAME = "frozen_input_manifest_v10_0_5_17.json"
CATALOG_NAME = "frozen_parameter_catalog_v10_0_5_17.json"
EXPECTED_KERNEL_FAMILY_SHA256 = (
    "a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde"
)
KERNEL_FAMILY_DESTINATION = (
    "signed_kernel/v10_2_14_active_only_campaign_family.json"
)

ENTRY_FILES = {
    "v10.2.25": (
        "arrhenius_fracture/sharp_front_v10_2_25_audited.py",
        "arrhenius_fracture/sharp_front_v10_2_25.py",
        "arrhenius_fracture/data/materials/v10_2_25_v913_paper_campaign_registry.csv",
        "arrhenius_fracture/data/materials/v10_2_25_v913_paper_campaign_selection.json",
    ),
    "v10.2.26": (
        "arrhenius_fracture/sharp_front_v10_2_26_audited.py",
        "arrhenius_fracture/sharp_front_v10_2_26.py",
        "arrhenius_fracture/data/materials/v10_2_26_v913_weakT_ceramic_registry.csv",
        "arrhenius_fracture/data/materials/v10_2_26_v913_weakT_ceramic_selection.json",
    ),
}

OPTIONS = (
    ("v10.2.25", "v913_paper_peak01_0242980_persistent_sites"),
    ("v10.2.25", "v913_paper_peak02_0127508_persistent_sites"),
    ("v10.2.25", "v913_paper_peak03_0115460_persistent_sites"),
    ("v10.2.25", "v913_paper_dbtt01_0202500_persistent_sites"),
    ("v10.2.25", "v913_paper_dbtt02_0088403_persistent_sites"),
    ("v10.2.26", "v913_paper_weakT01_0257068_persistent_sites"),
    ("v10.2.26", "v913_paper_ceramic01_0189364_persistent_sites"),
)


def _parse_configuration(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    result: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    required = {"release", "PF_repo_root", "PF_commit", "kernel_family"}
    missing = sorted(required.difference(result))
    if missing:
        raise ValueError(f"campaign configuration is missing fields {missing}")
    if result["release"] != "10.0.5.17":
        raise ValueError(f"unexpected campaign release {result['release']!r}")
    return result


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git_object_exists(repo: Path, expression: str) -> None:
    completed = subprocess.run(
        ["git", "-C", str(repo), "cat-file", "-e", expression],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"git object is unavailable in {repo}: {expression}\n"
            + completed.stderr.decode(errors="replace")
        )


def _git_file(repo: Path, commit: str, relative_path: str) -> bytes:
    completed = subprocess.run(
        ["git", "-C", str(repo), "show", f"{commit}:{relative_path}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"unable to read {relative_path!r} from commit {commit}\n"
            + completed.stderr.decode(errors="replace")
        )
    return completed.stdout


def _git_root(path: Path) -> Path | None:
    candidate = path.resolve() if path.is_dir() else path.resolve().parent
    completed = subprocess.run(
        ["git", "-C", str(candidate), "rev-parse", "--show-toplevel"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return None
    return Path(completed.stdout.strip()).resolve()


def _git_head(repo: Path | None) -> str | None:
    if repo is None:
        return None
    completed = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def _write_file(root: Path, relative_path: str, data: bytes) -> dict[str, Any]:
    destination = root / relative_path
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)
    return {
        "relative_path": relative_path,
        "sha256": _sha256_bytes(data),
        "bytes": len(data),
    }


def _validate_catalog(vendored_root: Path) -> list[dict[str, Any]]:
    records = []
    for entry, option in OPTIONS:
        candidate, audit = load_audited_parameter_option(vendored_root, entry, option)
        records.append(
            {
                "parameter_entry": entry,
                "parameter_option": option,
                "candidate_id": candidate.candidate_id,
                "material_class": audit.get("material_class"),
                "paper_role": audit.get("paper_role"),
                "selected_row_sha256": audit["selected_row_sha256"],
                "registry_sha256": audit["registry_sha256"],
                "selection_sha256": audit["selection_sha256"],
                "persistent_sites": audit["persistent_sites"],
                "finite_source_inventory": audit["finite_source_inventory"],
                "source_refresh": audit["source_refresh"],
                "explicit_recovery": audit["explicit_recovery"],
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument(
        "--kernel-family-source",
        type=Path,
        default=None,
        help="Override the kernel-family path recorded by the campaign.",
    )
    parser.add_argument(
        "--expected-kernel-family-sha256",
        default=EXPECTED_KERNEL_FAMILY_SHA256,
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    campaign_root = args.campaign_root.expanduser().resolve()
    configuration_path = campaign_root / CONFIGURATION_NAME
    configuration = _parse_configuration(configuration_path)

    parameter_source_repo = Path(configuration["PF_repo_root"]).expanduser().resolve()
    parameter_source_commit = configuration["PF_commit"]
    if not (parameter_source_repo / ".git").exists():
        raise FileNotFoundError(
            "recorded parameter-source repository is unavailable: "
            f"{parameter_source_repo}"
        )
    _git_object_exists(
        parameter_source_repo,
        f"{parameter_source_commit}^{{commit}}",
    )

    recorded_family_path = Path(configuration["kernel_family"]).expanduser().resolve()
    family_source = (
        args.kernel_family_source.expanduser().resolve()
        if args.kernel_family_source is not None
        else recorded_family_path
    )
    if not family_source.is_file():
        raise FileNotFoundError(
            "recorded kernel-family file is unavailable: "
            f"{family_source}; pass --kernel-family-source PATH only when using an "
            "identical byte-for-byte copy"
        )
    family_data = family_source.read_bytes()
    family_sha256 = _sha256_bytes(family_data)
    expected_family_sha256 = str(args.expected_kernel_family_sha256).strip()
    if expected_family_sha256 and family_sha256 != expected_family_sha256:
        raise ValueError(
            "kernel-family SHA mismatch: "
            f"{family_sha256} != {expected_family_sha256}"
        )

    required_parameter_paths = sorted(
        {path for paths in ENTRY_FILES.values() for path in paths}
    )

    destination = args.destination.expanduser().resolve()
    if destination.exists() and not args.force:
        raise FileExistsError(
            f"refusing to overwrite {destination}; pass --force for a clean replacement"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_parent = Path(
        tempfile.mkdtemp(prefix=destination.name + ".tmp.", dir=destination.parent)
    )
    staging = temporary_parent / destination.name
    staging.mkdir()

    try:
        file_records = []
        for relative_path in required_parameter_paths:
            data = _git_file(
                parameter_source_repo,
                parameter_source_commit,
                relative_path,
            )
            record = _write_file(staging, relative_path, data)
            record["source_kind"] = "parameter_source_commit"
            file_records.append(record)

        family_record = _write_file(
            staging,
            KERNEL_FAMILY_DESTINATION,
            family_data,
        )
        family_record["source_kind"] = "recorded_kernel_family"
        file_records.append(family_record)

        catalog_records = _validate_catalog(staging)
        family_repo = _git_root(family_source)
        manifest = {
            "schema": "v10.0.5.17_frozen_external_input_vendor_v2",
            "fem_czm_baseline_commit": BASELINE_COMMIT,
            "source_campaign_root": str(campaign_root),
            "source_campaign_configuration": str(configuration_path),
            "source_release": configuration["release"],
            "parameter_source_repository_recorded_path": str(parameter_source_repo),
            "parameter_source_repository_commit": parameter_source_commit,
            "parameter_files_vendored_from_recorded_commit": True,
            "kernel_family_recorded_path": str(recorded_family_path),
            "kernel_family_source_path": str(family_source),
            "kernel_family_source_repository": (
                str(family_repo) if family_repo is not None else None
            ),
            "kernel_family_source_repository_head": _git_head(family_repo),
            "kernel_family_relative_path": KERNEL_FAMILY_DESTINATION,
            "kernel_family_sha256": family_record["sha256"],
            "kernel_family_expected_sha256": expected_family_sha256,
            "kernel_family_vendored_byte_for_byte": True,
            "file_count": len(file_records),
            "files": file_records,
            "solver_or_constitutive_files_modified": False,
            "parameter_values_reconstructed": False,
        }
        (staging / MANIFEST_NAME).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        (staging / CATALOG_NAME).write_text(
            json.dumps(
                {
                    "schema": "v10.0.5.17_frozen_parameter_catalog_v1",
                    "parameter_source_repository_commit": parameter_source_commit,
                    "records": catalog_records,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

        if destination.exists():
            shutil.rmtree(destination)
        staging.replace(destination)
    finally:
        shutil.rmtree(temporary_parent, ignore_errors=True)

    result = {
        "destination": str(destination),
        "manifest": str(destination / MANIFEST_NAME),
        "catalog": str(destination / CATALOG_NAME),
        "parameter_source_commit": parameter_source_commit,
        "kernel_family": str(destination / KERNEL_FAMILY_DESTINATION),
        "kernel_family_sha256": family_record["sha256"],
        "validated_parameter_options": len(catalog_records),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
