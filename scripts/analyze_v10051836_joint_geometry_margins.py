#!/usr/bin/env python3
"""Combine v10.0.5.18.3.5/3.6 geometry audits into joint-margin tables."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable


DEFAULT_Q_FLOOR = 0.035
DEFAULT_A_FLOOR = 0.08


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def _classification(q: float | None, area: float | None, qfloor: float, afloor: float) -> str:
    if q is None or area is None:
        return "unavailable"
    qok = q >= qfloor
    aok = area >= afloor
    if qok and aok:
        return "joint_pass"
    if not qok and aok:
        return "quality_only_failure"
    if qok and not aok:
        return "area_only_failure"
    return "joint_failure"


def _row(
    *,
    source_file: Path,
    source_kind: str,
    event_index: int | None,
    refinement_level: int | None,
    partition_count: int | None,
    segment_index: int | None,
    q: Any,
    area: Any,
    qfloor: Any,
    afloor: Any,
    h_over_da: Any = None,
    clip_state: Any = None,
    minimum_barycentric: Any = None,
    distance_to_vertex_m: Any = None,
    distance_to_edge_m: Any = None,
    selected: Any = None,
    reason: Any = None,
) -> dict[str, Any]:
    qvalue = None if q is None else float(q)
    avalue = None if area is None else float(area)
    qf = float(DEFAULT_Q_FLOOR if qfloor is None else qfloor)
    af = float(DEFAULT_A_FLOOR if afloor is None else afloor)
    return {
        "source_file": str(source_file),
        "source_kind": source_kind,
        "event_index": event_index,
        "refinement_level": refinement_level,
        "partition_count": partition_count,
        "segment_index": segment_index,
        "min_triangle_quality": qvalue,
        "triangle_quality_floor": qf,
        "quality_margin": None if qvalue is None else qvalue / qf,
        "min_child_area_ratio": avalue,
        "child_area_ratio_floor": af,
        "area_margin": None if avalue is None else avalue / af,
        "limiting_margin": (
            None
            if qvalue is None or avalue is None
            else min(qvalue / qf, avalue / af)
        ),
        "joint_margin_classification": _classification(qvalue, avalue, qf, af),
        "committed_tip_h_over_da": h_over_da,
        "event_length_clip_state": clip_state,
        "minimum_barycentric_coordinate": minimum_barycentric,
        "distance_to_nearest_vertex_m": distance_to_vertex_m,
        "distance_to_nearest_edge_m": distance_to_edge_m,
        "selected": selected,
        "reason": reason,
    }


def _v36_rows(path: Path, payload: dict[str, Any]) -> Iterable[dict[str, Any]]:
    for index, event in enumerate(payload.get("events", [])):
        result = dict(event.get("backend_result", {}) or {})
        target = dict(event.get("target_parent", {}) or {})
        stochastic = dict(event.get("stochastic_event", {}) or {})
        yield _row(
            source_file=path,
            source_kind="v10051836_physical_event_backend_result",
            event_index=int(event.get("physical_event_index", index)),
            refinement_level=None,
            partition_count=None,
            segment_index=None,
            q=result.get("min_triangle_quality"),
            area=result.get("min_child_area_ratio"),
            qfloor=result.get("triangle_quality_floor"),
            afloor=result.get("child_area_ratio_floor"),
            h_over_da=event.get("committed_tip_h_over_requested_da"),
            clip_state=stochastic.get("event_length_clip_state"),
            minimum_barycentric=target.get("minimum_barycentric_coordinate"),
            distance_to_vertex_m=target.get("distance_to_nearest_vertex_m"),
            distance_to_edge_m=target.get("distance_to_nearest_edge_m"),
            selected=True,
            reason=result.get("reason"),
        )


def _v35_candidate_rows(
    path: Path,
    *,
    event_index: int | None,
    level: int | None,
    partition_count: int | None,
    segment_index: int | None,
    candidate: dict[str, Any],
    selected: bool,
    source_kind: str,
) -> Iterable[dict[str, Any]]:
    q = candidate.get(
        "predicted_min_triangle_quality",
        candidate.get("min_triangle_quality"),
    )
    area = candidate.get(
        "predicted_min_child_area_ratio",
        candidate.get("min_immediate_child_area_ratio"),
    )
    yield _row(
        source_file=path,
        source_kind=source_kind,
        event_index=event_index,
        refinement_level=level,
        partition_count=partition_count,
        segment_index=segment_index,
        q=q,
        area=area,
        qfloor=candidate.get("triangle_quality_floor"),
        afloor=candidate.get("child_area_ratio_floor"),
        selected=selected,
        reason=candidate.get("reason"),
    )


def _v35_rows(path: Path, payload: dict[str, Any]) -> Iterable[dict[str, Any]]:
    for index, veto in enumerate(payload.get("candidate_vetoes", [])):
        yield _row(
            source_file=path,
            source_kind="v10051835_candidate_veto",
            event_index=index,
            refinement_level=None,
            partition_count=None,
            segment_index=None,
            q=veto.get("min_triangle_quality"),
            area=veto.get("min_child_area_ratio"),
            qfloor=veto.get("triangle_quality_floor"),
            afloor=veto.get("child_area_ratio_floor"),
            h_over_da=veto.get("active_tip_h_over_da"),
            selected=False,
            reason=";".join(map(str, veto.get("issues", []))),
        )

    for index, refinement in enumerate(payload.get("patch_refinements", [])):
        level = refinement.get("level")
        partitions = refinement.get("physical_event_partition_count")
        segment = refinement.get("physical_event_segment_index")
        yield from _v35_candidate_rows(
            path,
            event_index=index,
            level=None if level is None else int(level),
            partition_count=None if partitions is None else int(partitions),
            segment_index=None if segment is None else int(segment),
            candidate=refinement,
            selected=bool(refinement.get("accepted", False)),
            source_kind="v10051835_selected_patch_refinement",
        )
        for error in refinement.get("errors", []) or []:
            if isinstance(error, dict):
                yield from _v35_candidate_rows(
                    path,
                    event_index=index,
                    level=None if level is None else int(level),
                    partition_count=None if partitions is None else int(partitions),
                    segment_index=None if segment is None else int(segment),
                    candidate=error,
                    selected=False,
                    source_kind="v10051835_rejected_patch_candidate",
                )


def _rows(path: Path, payload: dict[str, Any]) -> list[dict[str, Any]]:
    schema = str(payload.get("schema", ""))
    if "3.6" in schema or "committed_tip_resolution" in schema:
        return list(_v36_rows(path, payload))
    if "3.5" in schema or "quality_aware" in schema:
        return list(_v35_rows(path, payload))
    raise ValueError(f"unsupported audit schema in {path}: {schema!r}")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else ["source_file"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    classes: dict[str, int] = {}
    clip_states: dict[str, int] = {}
    for row in rows:
        classification = str(row["joint_margin_classification"])
        classes[classification] = classes.get(classification, 0) + 1
        clip = row.get("event_length_clip_state")
        if clip is not None:
            key = str(clip)
            clip_states[key] = clip_states.get(key, 0) + 1
    feasible = [row for row in rows if row["joint_margin_classification"] == "joint_pass"]
    return {
        "row_count": len(rows),
        "joint_margin_class_counts": classes,
        "event_length_clip_counts": clip_states,
        "joint_pass_row_count": len(feasible),
        "selected_joint_pass_row_count": sum(
            bool(row.get("selected")) for row in feasible
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("audits", nargs="+", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []
    for path in args.audits:
        rows.extend(_rows(path, _load(path)))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.out_dir / "joint_geometry_margins.csv", rows)
    payload = {
        "schema": "v10.0.5.18.3.6_joint_geometry_margin_analysis",
        "inputs": [str(path) for path in args.audits],
        "summary": _summary(rows),
        "rows": rows,
    }
    (args.out_dir / "joint_geometry_margins.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
    )
    print(json.dumps(payload["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
