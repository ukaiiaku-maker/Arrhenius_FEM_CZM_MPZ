#!/usr/bin/env python3
"""Run v9.12 candidates independently so one failure cannot abort a campaign."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


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
    return p.parse_args()


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

    print(f"RESILIENT_CAMPAIGN_START candidates={len(rows)} stage={a.stage}", flush=True)
    for i, row in enumerate(rows, start=1):
        cid = row.get("candidate_id", f"candidate_{i}")
        one = tmp / f"{cid}.csv"
        with one.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)
        cmd = [
            sys.executable, "-u", str(base),
            "--stage", a.stage,
            "--candidate-registry", str(one),
            "--protocol-csv", a.protocol_csv,
            "--physics-json", a.physics_json,
            "--temperatures", *a.temperatures,
            "--window-um", *a.window_um,
            "--target-cleavage-rate-s", a.target_cleavage_rate_s,
            "--min-amplitude", a.min_amplitude,
            "--target-localization", a.target_localization,
            "--max-width-K", a.max_width_K,
            "--out", str(out),
        ]
        print(f"RESILIENT_CANDIDATE_START index={i}/{len(rows)} candidate={cid}", flush=True)
        result = subprocess.run(cmd)
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
        records.append(rec)
        print(f"RESILIENT_CANDIDATE_RESULT candidate={cid} status={rec['status']}", flush=True)

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
    (out / "resilient_campaign_summary.json").write_text(
        json.dumps(
            {
                "candidate_count": len(records),
                "complete_count": sum(r.get("status") == "complete" for r in records),
                "unresolved_count": sum(r.get("status") != "complete" for r in records),
                "records": ordered,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"RESILIENT_CAMPAIGN_COMPLETE candidates={len(records)} out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
