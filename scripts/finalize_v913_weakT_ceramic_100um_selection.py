#!/usr/bin/env python3
"""Finalize the v9.13 100 um weak-T and ceramic paper selections.

The program reads the completed class-specific 1-D ranking tables, requires at
least two strict candidates in each class, retains the two lowest-score strict
rows, marks rank 1 as the paper primary and rank 2 as the backup, and exports a
hash-checked two-row primary handoff for the 2-D code.  Active candidate values
are copied exactly; no fitting, rounding, interpolation, or default substitution
is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)


CLASS_SPECS = {
    "weakT_FCC_like": {
        "ranked_file": "weakT_ranked.csv",
        "gate_column": "oneD_weakT_gate",
        "score_column": "oneD_weakT_score",
        "option_stem": "weakT",
    },
    "ceramic_like": {
        "ranked_file": "ceramic_ranked.csv",
        "gate_column": "oneD_ceramic_gate",
        "score_column": "oneD_ceramic_score",
        "option_stem": "ceramic",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--top-per-class", type=int, default=2)
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, bool):
        return value
    if hasattr(value, "item"):
        return json_safe(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n", ""}:
        return False
    raise ValueError(f"cannot interpret boolean value {value!r}")


def candidate_suffix(candidate_id: str) -> str:
    match = re.search(r"(\d+)$", candidate_id)
    if match is None:
        raise RuntimeError(f"candidate ID does not end in digits: {candidate_id}")
    return match.group(1)


def option_key(material_class: str, rank: int, candidate_id: str) -> str:
    stem = str(CLASS_SPECS[material_class]["option_stem"])
    suffix = candidate_suffix(candidate_id)
    return f"v913_paper_{stem}{rank:02d}_{suffix}_persistent_sites"


def validate_active_row(row: Mapping[str, Any]) -> dict[str, float]:
    candidate_id = str(row.get("candidate_id", ""))
    if not candidate_id:
        raise RuntimeError("candidate row lacks candidate_id")
    active = effective_candidate_parameters(row)
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        value = float(active[field])
        if not math.isfinite(value):
            raise RuntimeError(f"candidate {candidate_id} has nonfinite {field}={value}")
    return active


def select_class(
    frame: pd.DataFrame,
    material_class: str,
    count: int,
) -> pd.DataFrame:
    spec = CLASS_SPECS[material_class]
    required = {
        "candidate_id",
        "target_class",
        str(spec["gate_column"]),
        str(spec["score_column"]),
        *ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(
            f"{spec['ranked_file']} lacks required columns: {missing}"
        )
    local = frame[frame["target_class"].astype(str) == material_class].copy()
    if local.empty:
        raise RuntimeError(f"ranking table contains no {material_class} rows")
    local["_strict"] = local[str(spec["gate_column"])].map(as_bool)
    local["_score"] = pd.to_numeric(local[str(spec["score_column"])], errors="coerce")
    strict = local[local["_strict"] & local["_score"].map(math.isfinite)].copy()
    strict = strict.sort_values(
        ["_score", "candidate_id"], ascending=[True, True], kind="stable"
    ).reset_index(drop=True)
    if len(strict) < count:
        raise RuntimeError(
            f"{material_class} has only {len(strict)} strict candidates; "
            f"at least {count} are required"
        )
    selected = strict.head(count).copy()
    selected["paper_material_class"] = material_class
    selected["final_class_rank"] = range(1, len(selected) + 1)
    selected["selection_role"] = [
        "paper primary" if rank == 1 else "paper backup"
        for rank in selected["final_class_rank"]
    ]
    selected["option_key"] = [
        option_key(material_class, int(rank), str(candidate_id))
        for rank, candidate_id in zip(
            selected["final_class_rank"], selected["candidate_id"], strict=True
        )
    ]
    return selected


def metric_payload(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: json_safe(value)
        for key, value in row.items()
        if key.startswith("oneD_")
        and key not in {"oneD_weakT_gate", "oneD_ceramic_gate"}
    }


def build_outputs(
    analysis_dir: Path,
    top_per_class: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if top_per_class != 2:
        raise ValueError("the paper finalizer requires exactly two candidates per class")
    selections: list[pd.DataFrame] = []
    source_hashes: dict[str, str] = {}
    for material_class, spec in CLASS_SPECS.items():
        path = analysis_dir / str(spec["ranked_file"])
        if not path.is_file():
            raise FileNotFoundError(path)
        source_hashes[path.name] = sha256_path(path)
        frame = pd.read_csv(path)
        selections.append(select_class(frame, material_class, top_per_class))

    selected = pd.concat(selections, ignore_index=True, sort=False)
    records: list[dict[str, Any]] = []
    metadata: list[dict[str, Any]] = []
    for source in selected.to_dict(orient="records"):
        material_class = str(source["paper_material_class"])
        rank = int(source["final_class_rank"])
        active = validate_active_row(source)
        score_column = str(CLASS_SPECS[material_class]["score_column"])
        record = {
            "option_key": str(source["option_key"]),
            "candidate_id": str(source["candidate_id"]),
            "paper_material_class": material_class,
            "final_class_rank": rank,
            "selection_role": str(source["selection_role"]),
            "oneD_strict_gate_passed": True,
            "oneD_selection_score": float(source[score_column]),
            **active,
        }
        records.append(record)
        metadata.append(
            {
                "option_key": record["option_key"],
                "candidate_id": record["candidate_id"],
                "paper_material_class": material_class,
                "final_class_rank": rank,
                "selection_role": record["selection_role"],
                "oneD_strict_gate_passed": True,
                "oneD_selection_score": record["oneD_selection_score"],
                "oneD_metrics": metric_payload(source),
            }
        )

    columns = [
        "option_key",
        "candidate_id",
        "paper_material_class",
        "final_class_rank",
        "selection_role",
        "oneD_strict_gate_passed",
        "oneD_selection_score",
        *ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    ]
    top2 = pd.DataFrame(records)[columns]
    primary = top2[top2["final_class_rank"] == 1].copy().reset_index(drop=True)
    if len(primary) != 2:
        raise RuntimeError("expected exactly one primary candidate per class")
    manifest = {
        "schema": "v9.13_weakT_ceramic_100um_final_selection_v1",
        "source_analysis_dir": str(analysis_dir),
        "source_file_sha256": source_hashes,
        "target_extension_um": 100.0,
        "selection_rule": (
            "For each class, retain the two lowest-score rows that pass the strict "
            "100 um one-dimensional gate; rank 1 is the 2-D paper primary and rank 2 "
            "is the archived backup."
        ),
        "top_two_per_class": metadata,
        "primary_candidates": [row for row in metadata if row["final_class_rank"] == 1],
        "backup_candidates": [row for row in metadata if row["final_class_rank"] == 2],
        "active_parameter_fields": list(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
        "fixed_closure": {
            "persistent_sites": True,
            "finite_source_inventory": False,
            "source_depletion_on_emission": False,
            "source_refresh_on_crack_advance": False,
            "explicit_recovery": False,
            "dynamic_tip_radius": True,
            "dynamic_front_width": True,
        },
        "transfer_policy": (
            "Exact active-row transfer only; no fitting, transformation, rounding, "
            "or substitution of inactive legacy source/recovery coordinates."
        ),
    }
    return top2, primary, manifest


def main() -> int:
    args = parse_args()
    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.is_dir():
        raise FileNotFoundError(analysis_dir)
    top2, primary, manifest = build_outputs(analysis_dir, int(args.top_per_class))

    out = args.out_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    top2_path = out / "v9_13_weakT_ceramic_top2_registry.csv"
    primary_path = out / "v9_13_weakT_ceramic_primary_handoff.csv"
    manifest_path = out / "v9_13_weakT_ceramic_final_selection.json"
    top2.to_csv(top2_path, index=False)
    primary.to_csv(primary_path, index=False)

    top2_records = top2.to_dict(orient="records")
    primary_records = primary.to_dict(orient="records")
    manifest.update(
        {
            "top2_registry": str(top2_path),
            "top2_registry_sha256": sha256_path(top2_path),
            "top2_active_parameter_fingerprint_sha256": candidate_parameter_fingerprint(
                top2_records
            ),
            "primary_handoff_csv": str(primary_path),
            "primary_handoff_csv_sha256": sha256_path(primary_path),
            "primary_active_parameter_fingerprint_sha256": candidate_parameter_fingerprint(
                primary_records
            ),
        }
    )
    manifest_path.write_text(
        json.dumps(json_safe(manifest), indent=2, sort_keys=True) + "\n"
    )
    print("V913_WEAKT_CERAMIC_FINAL_SELECTION")
    for row in top2.to_dict(orient="records"):
        print(
            f"class={row['paper_material_class']} rank={row['final_class_rank']} "
            f"role={row['selection_role'].replace(' ', '_')} "
            f"candidate={row['candidate_id']} score={float(row['oneD_selection_score']):.9g}"
        )
    print(
        "V913_WEAKT_CERAMIC_PRIMARY_HANDOFF_WRITTEN "
        f"rows={len(primary)} csv={primary_path} manifest={manifest_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
