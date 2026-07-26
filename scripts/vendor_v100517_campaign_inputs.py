#!/usr/bin/env python3
"""Vendor the exact frozen inputs used by the v10.0.5.17 paper7 campaign.

This utility does not modify any FEM, CZM, stochastic-event, or moving-process-zone
code.  It reads the provenance written by the completed campaign and copies the
required parameter registries, stable option maps, audited-entry source records,
and signed-kernel family from the exact recorded source commit into this checkout.
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


def _relative_to_repo(path: Path, repo: Path) -> str:
    resolved = path.expanduser().resolve()
    try:
        return resolved.relative_to(repo.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError(
            "the recorded kernel family is not inside the recorded source repository: "
            f"family={resolved}, repository={repo.resolve()}"
        ) from exc


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
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    campaign_root = args.campaign_root.expanduser().resolve()
    configuration_path = campaign_root / CONFIGURATION_NAME
    configuration = _parse_configuration(configuration_path)

    source_repo = Path(configuration["PF_repo_root"]).expanduser().resolve()
    source_commit = configuration["PF_commit"]
    if not (source_repo / ".git").exists():
        raise FileNotFoundError(f"recorded source repository is unavailable: {source_repo}")
    _git_object_exists(source_repo, f"{source_commit}^{{commit}}")

    family_relative = _relative_to_repo(Path(configuration["kernel_family"]), source_repo)
    required_paths = sorted(
        {path for paths in ENTRY_FILES.values() for path in paths} | {family_relative}
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
        for relative_path in required_paths:
            data = _git_file(source_repo, source_commit, relative_path)
            file_records.append(_write_file(staging, relative_path, data))

        catalog_records = _validate_catalog(staging)
        family_record = next(
            row for row in file_records if row["relative_path"] == family_relative
        )
        manifest = {
            "schema": "v10.0.5.17_frozen_external_input_vendor_v1",
            "fem_czm_baseline_commit": BASELINE_COMMIT,
            "source_campaign_root": str(campaign_root),
            "source_campaign_configuration": str(configuration_path),
            "source_release": configuration["release"],
            "source_repository_recorded_path": str(source_repo),
            "source_repository_commit": source_commit,
            "kernel_family_relative_path": family_relative,
            "kernel_family_sha256": family_record["sha256"],
            "file_count": len(file_records),
            "files": file_records,
            "solver_or_constitutive_files_modified": False,
            "parameter_values_reconstructed": False,
            "vendored_byte_for_byte_from_recorded_commit": True,
        }
        (staging / MANIFEST_NAME).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        (staging / CATALOG_NAME).write_text(
            json.dumps(
                {
                    "schema": "v10.0.5.17_frozen_parameter_catalog_v1",
                    "source_repository_commit": source_commit,
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
        "source_commit": source_commit,
        "kernel_family": str(destination / family_relative),
        "kernel_family_sha256": family_record["sha256"],
        "validated_parameter_options": len(catalog_records),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
