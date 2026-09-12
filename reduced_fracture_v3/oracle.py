"""Read-only normalization of retained V5 evidence into a V3 oracle table."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping


ORACLE_SCHEMA = "oneD.v3.v5-oracle-export/1"
PINNED_V5_REPOSITORY = "ukaiiaku-maker/PF-fracture-fatigue"
PINNED_V5_HARDENING_SHA = "b58997bdb18cf4e9a32c251c073115d8b405bb27"
PINNED_V5_ATTESTATION_SHA = "c7583ecd0a259f28ce92780833d21a358a920f45"

MAP_REQUIRED_FIELDS = (
    "G_kinetic_used_J_per_m2",
    "cavity_sigma_nn_Pa",
    "cavity_sigma_tt_Pa",
    "cavity_tau_nt_Pa",
    "cavity_sigma_m_Pa",
    "compliance_m2_per_N",
    "potential_energy_J_per_m",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _repository_path(path: Path) -> str:
    parts = path.parts
    if "artifacts" in parts:
        return Path(*parts[parts.index("artifacts"):]).as_posix()
    return path.name


def _checkpoint_phase(void_state: Mapping[str, Any]) -> str:
    cavities = void_state.get("cavities", [])
    if cavities:
        return str(cavities[0]["phase"])
    sites = void_state.get("sites", [])
    return "NO_VOID" if not sites else str(sites[0]["phase"])


def _active_tip_position(checkpoint: Mapping[str, Any]) -> float | None:
    active = set(checkpoint.get("active_tip_ids", []))
    for branch in checkpoint.get("crack_network", {}).get("branches", []):
        if branch.get("branch_id") in active or branch.get("status") == "active":
            path = branch.get("path_m", [])
            if path:
                return float(path[-1][0])
    return None


def _base_row(kind: str, source_path: Path, source_sha: str) -> dict[str, Any]:
    return {
        "record_kind": kind,
        "source_path": _repository_path(source_path),
        "source_file_sha256": _sha256(source_path),
        "source_implementation_sha": source_sha,
        "oracle_repository": PINNED_V5_REPOSITORY,
        "oracle_hardening_sha": PINNED_V5_HARDENING_SHA,
        "oracle_attestation_sha": PINNED_V5_ATTESTATION_SHA,
    }


def static_rows(path: str | Path) -> list[dict[str, Any]]:
    source = Path(path)
    payload = json.loads(source.read_text())
    source_sha = str(payload["implementation_git_sha"])
    rows = []
    for item in payload["rows"]:
        configuration = item.get("configuration", {})
        observables = item.get("observables", {})
        crack_path = configuration.get("crack_path_m_requested") or []
        tip = crack_path[-1][0] if crack_path else None
        center = configuration.get("cavity_center_m") or [None, None]
        radius = configuration.get("cavity_radius_m")
        ligament = None
        if tip is not None and center[0] is not None and radius is not None:
            ligament = float(center[0]) - float(radius) - float(tip)
        row = {
            **_base_row("STATIC_MECHANICS", source, source_sha),
            "record_id": f"static:{item['case']}",
            "accepted_state_identity": None,
            "phase": "PRECONNECTION" if configuration.get("crack_enabled") else "VOID_ONLY",
            "x_tip_m": tip,
            "x_void_m": center[0],
            "lateral_offset_m": center[1],
            "radius_m": radius,
            "ligament_m": ligament,
            "load": configuration.get("opening_m"),
            "candidate_id": None,
            "G_local_J_per_m2": None,
            "G_marginal_J_per_m2": None,
            "G_kinetic_used_J_per_m2": None,
            "K_energy_equivalent_Pa_sqrt_m": None,
            "cavity_sigma_nn_Pa": None,
            "cavity_sigma_tt_Pa": None,
            "cavity_tau_nt_Pa": None,
            "cavity_sigma_m_Pa": None,
            "reaction_N_per_m": observables.get("reaction_top_N_per_m"),
            "compliance_m2_per_N": observables.get("compliance_m2_per_N"),
            "potential_energy_J_per_m": observables.get("stored_energy_J_per_m"),
            "topology_identity": item.get("actual_geometry_fingerprint"),
            "length_ledgers": None,
        }
        row["missing_map_fields"] = [name for name in MAP_REQUIRED_FIELDS if row.get(name) is None]
        row["map_qualification"] = "MAP_READY" if not row["missing_map_fields"] else "INCOMPLETE_RETAINED_EVIDENCE"
        rows.append(row)
    return rows


def checkpoint_rows(directory: str | Path) -> list[dict[str, Any]]:
    root = Path(directory)
    campaign_manifest = root.parent / "campaign_manifest.json"
    checkpoint_source_sha = PINNED_V5_ATTESTATION_SHA
    if campaign_manifest.exists():
        checkpoint_source_sha = str(json.loads(campaign_manifest.read_text())["implementation_sha"])
    rows = []
    for source in sorted(root.glob("*.json")):
        checkpoint = json.loads(source.read_text())
        void = checkpoint.get("production_void_state") or {}
        cavities = void.get("cavities", [])
        sites = void.get("sites", [])
        cavity = cavities[0] if cavities else None
        site = sites[0] if sites else None
        center = (cavity or site or {}).get("center_m", [None, None])
        radius = None if cavity is None else cavity.get("radius_m")
        tip = _active_tip_position(checkpoint)
        energy = checkpoint.get("energy_ledgers", {})
        row = {
            **_base_row("ACCEPTED_CHECKPOINT", source, checkpoint_source_sha),
            "record_id": f"checkpoint:{source.stem}",
            "accepted_state_identity": checkpoint.get("state_sha256"),
            "phase": _checkpoint_phase(void),
            "x_tip_m": tip,
            "x_void_m": center[0],
            "lateral_offset_m": center[1],
            "radius_m": radius,
            "ligament_m": (
                None if tip is None or radius is None else float(center[0]) - float(radius) - tip
            ),
            "load": None,
            "candidate_id": None,
            "G_local_J_per_m2": None,
            "G_marginal_J_per_m2": None,
            "G_kinetic_used_J_per_m2": None,
            "K_energy_equivalent_Pa_sqrt_m": None,
            "cavity_sigma_nn_Pa": None,
            "cavity_sigma_tt_Pa": None,
            "cavity_tau_nt_Pa": None,
            "cavity_sigma_m_Pa": None,
            "reaction_N_per_m": energy.get("latest_energy_reaction_identity"),
            "compliance_m2_per_N": None,
            "potential_energy_J_per_m": energy.get("latest_fem_energy_J_per_m"),
            "topology_identity": checkpoint.get("v12_support_state", {}).get(
                "complete_crack_graph_fingerprint"
            ),
            "length_ledgers": void.get("length_ledgers"),
            "source_v5_void_state_sha256": checkpoint.get("production_void_state_sha256"),
            "source_topology_transaction_model_id": checkpoint.get(
                "topology_transaction_model_id"
            ),
        }
        row["missing_map_fields"] = [name for name in MAP_REQUIRED_FIELDS if row.get(name) is None]
        row["map_qualification"] = "MAP_READY" if not row["missing_map_fields"] else "STATE_TRANSFER_ONLY"
        rows.append(row)
    return rows


def readiness(rows: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    rows = tuple(rows)
    missing = {name: sum(row.get(name) is None for row in rows) for name in MAP_REQUIRED_FIELDS}
    map_rows = tuple(row for row in rows if row.get("map_qualification") == "MAP_READY")
    preconnection = tuple(row for row in map_rows if row.get("phase") == "PRECONNECTION")
    radii = {float(row["radius_m"]) for row in preconnection}
    ligament_ratios = {
        float(row["ligament_m"]) / float(row["radius_m"])
        for row in preconnection
        if row.get("ligament_m") is not None and float(row["radius_m"]) > 0.0
    }
    geometry_points = {
        (float(row["ligament_m"]) / float(row["radius_m"]), float(row["radius_m"]))
        for row in preconnection
        if row.get("ligament_m") is not None and float(row["radius_m"]) > 0.0
    }
    preconnection_grid_complete = bool(
        len(radii) >= 3
        and len(ligament_ratios) >= 2
        and len(geometry_points) == len(radii) * len(ligament_ratios)
    )
    other_regimes = {str(row.get("phase")) for row in map_rows}
    topology_coverage_complete = {
        "PRECONNECTION",
        "CONNECTED_VOID",
        "DOWNSTREAM_FRONT_ACTIVE",
    }.issubset(other_regimes)
    return {
        "total_records": len(rows),
        "map_ready_records": len(map_rows),
        "missing_field_counts": missing,
        "preconnection_complete_geometry_grid": preconnection_grid_complete,
        "all_three_topology_regimes_present": topology_coverage_complete,
        "classification": (
            "READY_FOR_HELD_OUT_FIT"
            if preconnection_grid_complete and topology_coverage_complete
            else "BLOCKED_ORACLE_COORDINATES_INCOMPLETE"
        ),
    }


def export_v5_oracle(
    *,
    static_evidence: str | Path,
    checkpoint_directory: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    rows = static_rows(static_evidence) + checkpoint_rows(checkpoint_directory)
    payload = {
        "schema": ORACLE_SCHEMA,
        "source_policy": "PINNED_READ_ONLY_RETAINED_EVIDENCE_NO_NEW_2D_SOLVE",
        "pinned_v5_repository": PINNED_V5_REPOSITORY,
        "pinned_v5_hardening_sha": PINNED_V5_HARDENING_SHA,
        "pinned_v5_attestation_sha": PINNED_V5_ATTESTATION_SHA,
        "readiness": readiness(rows),
        "rows": rows,
    }
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


__all__ = [
    "MAP_REQUIRED_FIELDS",
    "ORACLE_SCHEMA",
    "PINNED_V5_ATTESTATION_SHA",
    "PINNED_V5_HARDENING_SHA",
    "checkpoint_rows",
    "export_v5_oracle",
    "readiness",
    "static_rows",
]
