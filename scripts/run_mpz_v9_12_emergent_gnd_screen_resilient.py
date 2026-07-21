#!/usr/bin/env python3
"""Run v9.12 candidates independently so one failure cannot abort a campaign."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time


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


def main() -> int:
    a = args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with Path(a.candidate_registry).open(newline="") as fp:
        reader = csv.DictReader(fp)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    tmp = out / "_single_candidate_registries"
    tmp.mkdir(exist_ok=True)
    records = []
    base = Path(__file__).with_name("run_mpz_v9_12_emergent_gnd_screen.py")
    progress_path = out / "resilient_progress.json"
    started_at_utc = utc_now()
    campaign_start = time.perf_counter()
    total = len(rows)
    progress_every = max(int(a.progress_every), 1)

    print(
        "RESILIENT_CAMPAIGN_START "
        f"candidates={total} stage={a.stage} started_at={started_at_utc}",
        flush=True,
    )
    write_progress(
        progress_path,
        progress_payload(
            started_at_utc=started_at_utc,
            completed=0,
            total=total,
            complete_count=0,
            unresolved_count=0,
            elapsed_s=0.0,
        ),
    )

    for i, row in enumerate(rows, start=1):
        cid = row.get("candidate_id", f"candidate_{i}")
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
        print(
            f"RESILIENT_CANDIDATE_START index={i}/{total} candidate={cid}",
            flush=True,
        )
        candidate_start = time.perf_counter()
        result = subprocess.run(cmd)
        candidate_elapsed_s = time.perf_counter() - candidate_start
        ranking = out / "ranking.csv"
        if result.returncode == 0 and ranking.exists():
            with ranking.open(newline="") as fp:
                rec = next(csv.DictReader(fp))
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
        rec["candidate_elapsed_s"] = f"{candidate_elapsed_s:.9g}"
        records.append(rec)
        print(
            "RESILIENT_CANDIDATE_RESULT "
            f"candidate={cid} status={rec['status']} "
            f"candidate_elapsed_s={candidate_elapsed_s:.3f}",
            flush=True,
        )

        elapsed_s = time.perf_counter() - campaign_start
        complete_count = sum(r.get("status") == "complete" for r in records)
        unresolved_count = len(records) - complete_count
        payload = progress_payload(
            started_at_utc=started_at_utc,
            completed=i,
            total=total,
            complete_count=complete_count,
            unresolved_count=unresolved_count,
            elapsed_s=elapsed_s,
        )
        write_progress(progress_path, payload)
        if i % progress_every == 0 or i == total:
            eta = payload["eta_s"]
            eta_text = "unknown" if eta is None else f"{float(eta):.1f}"
            print(
                "RESILIENT_PROGRESS "
                f"completed={i}/{total} complete={complete_count} "
                f"unresolved={unresolved_count} elapsed_s={elapsed_s:.1f} "
                f"mean_candidate_s={float(payload['mean_candidate_s']):.3f} "
                f"candidates_per_min={float(payload['candidates_per_min']):.3f} "
                f"eta_s={eta_text}",
                flush=True,
            )

    all_fields = ["rank"]
    for rec in records:
        for key in rec:
            if key not in all_fields:
                all_fields.append(key)
    ordered = sorted(
        records,
        key=lambda r: (
            r.get("status") == "complete",
            float(r.get("score", "-1e300")),
        ),
        reverse=True,
    )
    with (out / "ranking.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=all_fields)
        writer.writeheader()
        for rank, rec in enumerate(ordered, start=1):
            writer.writerow({"rank": rank, **rec})

    total_elapsed_s = time.perf_counter() - campaign_start
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
                "records": ordered,
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
