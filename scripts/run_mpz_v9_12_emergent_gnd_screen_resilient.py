#!/usr/bin/env python3
"""Run v9.12 candidates independently so one failure cannot abort a campaign."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time


CHECKPOINT_FIELDS = [
    "candidate_id",
    "stage",
    "status",
    "failure_reason",
    "candidate_elapsed_s",
    "score",
    "pass",
    "amplitude_MPa_sqrt_m",
    "largest_jump_localization",
    "transition_width_10_90_K",
    "linear_r2",
    "max_abs_K_shield_MPa_sqrt_m",
    "max_tau_gnd_tip_MPa",
    "max_gnd_abs_line_count_per_unit_thickness",
    "min_source_available_fraction",
]


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-registry", required=True)
    p.add_argument("--protocol-csv", required=True)
    p.add_argument("--physics-json", required=True)
    p.add_argument("--temperatures", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--stage", choices=("0d", "1d"), default="1d")
    p.add_argument("--window-um", nargs=2, default=("10", "30"))
    p.add_argument("--target-cleavage-rate-s", default="1e-3")
    p.add_argument("--min-amplitude", default="8")
    p.add_argument("--target-localization", default="0.5")
    p.add_argument("--max-width-K", default="200")
    p.add_argument(
        "--progress-every",
        type=int,
        default=1,
        help="Print aggregate progress every N completed candidates.",
    )
    p.add_argument(
        "--compact-output",
        action="store_true",
        help="Omit per-temperature JSON files while retaining candidate summaries.",
    )
    p.add_argument(
        "--quiet-inner",
        action="store_true",
        help="Suppress successful child-screen output; retain it only for failures.",
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume from resilient_records_checkpoint.csv in the output directory.",
    )
    p.add_argument(
        "--keep-single-registries",
        action="store_true",
        help="Keep temporary one-row registries after each candidate.",
    )
    return p.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def progress_payload(
    *,
    started_at_utc: str,
    completed: int,
    total: int,
    complete_count: int,
    unresolved_count: int,
    elapsed_s: float,
) -> dict[str, object]:
    rate_per_s = completed / elapsed_s if elapsed_s > 0.0 else 0.0
    eta_s = (total - completed) / rate_per_s if rate_per_s > 0.0 else None
    return {
        "started_at_utc": started_at_utc,
        "updated_at_utc": utc_now(),
        "completed_candidates": completed,
        "total_candidates": total,
        "complete_count": complete_count,
        "unresolved_count": unresolved_count,
        "elapsed_s": elapsed_s,
        "mean_candidate_s": elapsed_s / completed if completed else None,
        "candidates_per_min": 60.0 * rate_per_s,
        "eta_s": eta_s,
    }


def write_progress(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def read_csv_records(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as fp:
        return list(csv.DictReader(fp))


def append_checkpoint(path: Path, record: dict[str, object]) -> None:
    new_file = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=CHECKPOINT_FIELDS,
            extrasaction="ignore",
        )
        if new_file:
            writer.writeheader()
        writer.writerow(record)
        fp.flush()
        os.fsync(fp.fileno())


def record_fields(records: list[dict[str, object]], *, include_rank: bool) -> list[str]:
    fields = ["rank"] if include_rank else []
    for key in CHECKPOINT_FIELDS:
        if any(key in rec for rec in records) and key not in fields:
            fields.append(key)
    for rec in records:
        for key in rec:
            if key != "rank" and key not in fields:
                fields.append(key)
    return fields


def ordered_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        records,
        key=lambda r: (
            r.get("status") == "complete",
            float(r.get("score", "-1e300")),
        ),
        reverse=True,
    )


def write_ranking(path: Path, records: list[dict[str, object]]) -> None:
    ordered = ordered_records(records)
    fields = record_fields(ordered, include_rank=True)
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for rank, rec in enumerate(ordered, start=1):
            clean = {key: value for key, value in rec.items() if key != "rank"}
            writer.writerow({"rank": rank, **clean})


def main() -> int:
    a = args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with Path(a.candidate_registry).open(newline="") as fp:
        reader = csv.DictReader(fp)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if not rows:
        raise RuntimeError(f"candidate registry is empty: {a.candidate_registry}")

    tmp = out / "_single_candidate_registries"
    tmp.mkdir(exist_ok=True)
    checkpoint_path = out / "resilient_records_checkpoint.csv"
    progress_path = out / "resilient_progress.json"
    base = Path(__file__).with_name("run_mpz_v9_12_emergent_gnd_screen.py")

    if a.resume:
        records: list[dict[str, object]] = list(read_csv_records(checkpoint_path))
        if not records:
            print(
                f"RESILIENT_RESUME_NO_CHECKPOINT path={checkpoint_path}; starting fresh",
                flush=True,
            )
    else:
        if checkpoint_path.exists():
            raise RuntimeError(
                f"checkpoint already exists: {checkpoint_path}; "
                "use --resume or choose a new output directory"
            )
        records = []

    done_ids = {str(rec.get("candidate_id", "")) for rec in records}
    pending = [
        row
        for row in rows
        if str(row.get("candidate_id", "")) not in done_ids
    ]
    prior_elapsed_s = sum(
        float(rec.get("candidate_elapsed_s", 0.0) or 0.0) for rec in records
    )
    started_at_utc = utc_now()
    campaign_start = time.perf_counter()
    total = len(rows)
    progress_every = max(int(a.progress_every), 1)

    print(
        "RESILIENT_CAMPAIGN_START "
        f"candidates={total} already_completed={len(records)} "
        f"pending={len(pending)} stage={a.stage} started_at={started_at_utc}",
        flush=True,
    )
    initial_complete = sum(r.get("status") == "complete" for r in records)
    write_progress(
        progress_path,
        progress_payload(
            started_at_utc=started_at_utc,
            completed=len(records),
            total=total,
            complete_count=initial_complete,
            unresolved_count=len(records) - initial_complete,
            elapsed_s=prior_elapsed_s,
        ),
    )

    for row in pending:
        completed_index = len(records) + 1
        cid = str(row.get("candidate_id", f"candidate_{completed_index}"))
        one = tmp / f"{cid}.csv"
        with one.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)

        cmd = [
            sys.executable,
            "-u",
            str(base),
            "--stage",
            a.stage,
            "--candidate-registry",
            str(one),
            "--protocol-csv",
            a.protocol_csv,
            "--physics-json",
            a.physics_json,
            "--temperatures",
            *a.temperatures,
            "--window-um",
            *a.window_um,
            "--target-cleavage-rate-s",
            a.target_cleavage_rate_s,
            "--min-amplitude",
            a.min_amplitude,
            "--target-localization",
            a.target_localization,
            "--max-width-K",
            a.max_width_K,
            "--out",
            str(out),
        ]
        if a.compact_output:
            cmd.append("--compact-output")
        if a.quiet_inner:
            cmd.append("--quiet-cases")

        print(
            f"RESILIENT_CANDIDATE_START index={completed_index}/{total} candidate={cid}",
            flush=True,
        )
        candidate_start = time.perf_counter()
        if a.quiet_inner:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
        else:
            result = subprocess.run(cmd)
        candidate_elapsed_s = time.perf_counter() - candidate_start

        child_ranking = out / "ranking.csv"
        if result.returncode == 0 and child_ranking.exists():
            with child_ranking.open(newline="") as fp:
                rec: dict[str, object] = dict(next(csv.DictReader(fp)))
            rec.pop("rank", None)
            rec["status"] = "complete"
            rec["failure_reason"] = ""
        else:
            rec = {
                "candidate_id": cid,
                "stage": a.stage,
                "status": "unresolved",
                "failure_reason": f"single-candidate returncode={result.returncode}",
                "score": "-1e300",
                "pass": "False",
            }
            if a.quiet_inner:
                failure_root = out / cid
                failure_root.mkdir(parents=True, exist_ok=True)
                (failure_root / "candidate_failure.log").write_text(
                    result.stdout or ""
                )
                tail = " | ".join((result.stdout or "").splitlines()[-3:])
                if tail:
                    print(
                        f"RESILIENT_FAILURE_TAIL candidate={cid} {tail}",
                        flush=True,
                    )

        rec["candidate_elapsed_s"] = f"{candidate_elapsed_s:.9g}"
        records.append(rec)
        append_checkpoint(checkpoint_path, rec)
        if not a.keep_single_registries:
            one.unlink(missing_ok=True)

        print(
            "RESILIENT_CANDIDATE_RESULT "
            f"candidate={cid} status={rec['status']} "
            f"candidate_elapsed_s={candidate_elapsed_s:.3f}",
            flush=True,
        )

        elapsed_s = prior_elapsed_s + (time.perf_counter() - campaign_start)
        complete_count = sum(r.get("status") == "complete" for r in records)
        unresolved_count = len(records) - complete_count
        payload = progress_payload(
            started_at_utc=started_at_utc,
            completed=len(records),
            total=total,
            complete_count=complete_count,
            unresolved_count=unresolved_count,
            elapsed_s=elapsed_s,
        )
        write_progress(progress_path, payload)
        if len(records) % progress_every == 0 or len(records) == total:
            eta = payload["eta_s"]
            eta_text = "unknown" if eta is None else f"{float(eta):.1f}"
            print(
                "RESILIENT_PROGRESS "
                f"completed={len(records)}/{total} complete={complete_count} "
                f"unresolved={unresolved_count} elapsed_s={elapsed_s:.1f} "
                f"mean_candidate_s={float(payload['mean_candidate_s']):.3f} "
                f"candidates_per_min={float(payload['candidates_per_min']):.3f} "
                f"eta_s={eta_text}",
                flush=True,
            )

    write_ranking(out / "ranking.csv", records)
    total_elapsed_s = prior_elapsed_s + (time.perf_counter() - campaign_start)
    complete_count = sum(r.get("status") == "complete" for r in records)
    unresolved_count = len(records) - complete_count
    final_progress = progress_payload(
        started_at_utc=started_at_utc,
        completed=len(records),
        total=total,
        complete_count=complete_count,
        unresolved_count=unresolved_count,
        elapsed_s=total_elapsed_s,
    )
    write_progress(progress_path, final_progress)
    (out / "resilient_campaign_summary.json").write_text(
        json.dumps(
            {
                **final_progress,
                "candidate_count": len(records),
                "checkpoint_csv": str(checkpoint_path.resolve()),
                "compact_output": bool(a.compact_output),
                "records": ordered_records(records),
            },
            indent=2,
        )
        + "\n"
    )
    print(
        "RESILIENT_CAMPAIGN_COMPLETE "
        f"candidates={len(records)} complete={complete_count} "
        f"unresolved={unresolved_count} elapsed_s={total_elapsed_s:.3f} "
        f"candidates_per_min={float(final_progress['candidates_per_min']):.3f} "
        f"out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
