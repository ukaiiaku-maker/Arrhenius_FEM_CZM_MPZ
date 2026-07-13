#!/usr/bin/env python3
"""Read-only live monitor for mixed-mode FEM/CZM v2 event-controlled campaigns.

This script uses only the Python standard library. It does not import NumPy,
SciPy, or arrhenius_fracture and never writes inside the campaign directory.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

STEP_RE = re.compile(
    r"\[T=(?P<T>[-+0-9.eE]+)K\]\s+step\s+(?P<step>\d+)\s+"
    r"KJ=\s*(?P<KJ>[-+0-9.eE]+)\s+"
    r"sig_tip=\s*(?P<sig>[-+0-9.eE]+)GPa\s+"
    r"B=\s*(?P<B>[-+0-9.eE]+)\s+"
    r"N_em=\s*(?P<Nem>[-+0-9.eE]+)\s+"
    r"a=(?P<a>[-+0-9.eE]+)mm\s+nfr=(?P<nfr>\d+)"
)
PSI_TAG_RE = re.compile(r"psi_([mp])(\d+)p(\d+)$")
ITER_RE = re.compile(r"iter_(\d+)_alpha_([mp])(\d+)p(\d+)$")


def parse_items(text: str, cast=str) -> list[Any]:
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def parse_psi_tag(name: str) -> float | None:
    m = PSI_TAG_RE.fullmatch(name)
    if not m:
        return None
    sign = -1.0 if m.group(1) == "m" else 1.0
    return sign * (int(m.group(2)) + int(m.group(3)) / 10.0)


def psi_tag(value: float) -> str:
    return "psi_" + ("m" if value < 0 else "p") + f"{abs(value):05.1f}".replace(".", "p")


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def read_last_lines(path: Path, max_bytes: int = 65536) -> list[str]:
    try:
        with path.open("rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - max_bytes), os.SEEK_SET)
            data = fh.read().decode("utf-8", errors="replace")
        return data.splitlines()
    except OSError:
        return []


def parse_max_steps(command_path: Path, fallback: int) -> int:
    try:
        tokens = shlex.split(command_path.read_text())
        idx = tokens.index("--steps")
        return int(tokens[idx + 1])
    except (OSError, ValueError, IndexError):
        return fallback


def latest_step(log_path: Path) -> dict[str, Any] | None:
    found = None
    for line in read_last_lines(log_path):
        m = STEP_RE.search(line)
        if m:
            found = {
                "step": int(m.group("step")),
                "KJ": float(m.group("KJ")),
                "sig_tip": float(m.group("sig")),
                "B": float(m.group("B")),
                "N_em": float(m.group("Nem")),
                "a_mm": float(m.group("a")),
                "nfr": int(m.group("nfr")),
                "advanced": "ADVANCE" in line,
            }
    return found


def fmt_num(value: Any, width: int = 7, decimals: int = 1) -> str:
    try:
        x = float(value)
        if not math.isfinite(x):
            return "-".rjust(width)
        return f"{x:{width}.{decimals}f}"
    except (TypeError, ValueError):
        return "-".rjust(width)


def fmt_bool(value: Any) -> str:
    if value is True:
        return "yes"
    if value is False:
        return "no"
    return "-"


def seconds_since(path: Path) -> float:
    try:
        return max(0.0, time.time() - path.stat().st_mtime)
    except OSError:
        return float("inf")


def age_text(seconds: float) -> str:
    if not math.isfinite(seconds):
        return "-"
    if seconds < 60:
        return f"{int(seconds)}s"
    if seconds < 3600:
        return f"{int(seconds // 60)}m"
    return f"{seconds / 3600:.1f}h"


@dataclass
class CaseStatus:
    klass: str
    target: float
    state: str
    iteration: int | None = None
    alpha: float | None = None
    achieved: float | None = None
    error: float | None = None
    reliable: bool | None = None
    mode: str = "-"
    step: int | None = None
    max_steps: int | None = None
    KJ: float | None = None
    sig_tip: float | None = None
    age_s: float = float("inf")
    trial_dir: str = ""
    note: str = ""


def infer_latest_trial(case_root: Path) -> Path | None:
    trials = [p for p in (case_root / "trials").glob("iter_*_alpha_*") if p.is_dir()]
    if not trials:
        return None
    def key(p: Path) -> tuple[int, float]:
        m = ITER_RE.fullmatch(p.name)
        it = int(m.group(1)) if m else -1
        try:
            mt = p.stat().st_mtime
        except OSError:
            mt = 0.0
        return it, mt
    return max(trials, key=key)


def status_for_case(root: Path, klass: str, target: float, max_steps_default: int, stale_s: float) -> CaseStatus:
    case_root = root / klass / psi_tag(target)
    final_json = case_root / "mixed_mode_v2_final_summary.json"
    final = load_json(final_json)
    if final is not None:
        converged = bool(final.get("event_state_control_converged", False))
        return CaseStatus(
            klass=klass,
            target=target,
            state="CONVERGED" if converged else "NOT_CONV",
            iteration=int(final.get("control_iterations", 0)) - 1,
            alpha=final.get("selected_loading_angle_deg"),
            achieved=final.get("achieved_psi_deg"),
            error=final.get("psi_error_deg"),
            reliable=final.get("projection_reliable"),
            mode=str(final.get("mode_classification", "-")),
            KJ=final.get("KJ_first_MPa_sqrt_m"),
            age_s=seconds_since(final_json),
            trial_dir=str(final.get("selected_trial_dir", "")),
        )

    trial = infer_latest_trial(case_root)
    if trial is None:
        return CaseStatus(klass, target, "PENDING")

    meta = load_json(trial / "control_trial_metadata.json") or {}
    summary = load_json(trial / "mixed_mode_first_passage_summary.json")
    log = trial / "run.log"
    step_data = latest_step(log)
    max_steps = parse_max_steps(trial / "command.txt", max_steps_default)
    log_age = seconds_since(log)
    m = ITER_RE.fullmatch(trial.name)
    iteration = int(meta.get("iteration", m.group(1) if m else 0))
    alpha = meta.get("loading_angle_deg")

    if summary is not None:
        achieved = summary.get("mode_phase_first_deg")
        err = None
        try:
            err = (float(achieved) - float(target) + 180.0) % 360.0 - 180.0
        except (TypeError, ValueError):
            pass
        return CaseStatus(
            klass, target, "TRIAL_DONE", iteration, alpha, achieved, err,
            summary.get("projection_reliable"), str(summary.get("mode_classification", "-")),
            step_data.get("step") if step_data else None, max_steps,
            summary.get("KJ_first_MPa_sqrt_m") or (step_data.get("KJ") if step_data else None),
            step_data.get("sig_tip") if step_data else None,
            log_age, str(trial), "controller update pending",
        )

    if not log.exists():
        state = "STARTING"
    elif log_age > stale_s:
        state = "STALE?"
    else:
        state = "RUNNING"
    return CaseStatus(
        klass, target, state, iteration, alpha,
        step=step_data.get("step") if step_data else None,
        max_steps=max_steps,
        KJ=step_data.get("KJ") if step_data else None,
        sig_tip=step_data.get("sig_tip") if step_data else None,
        age_s=log_age,
        trial_dir=str(trial),
        note="no recent log update" if state == "STALE?" else "",
    )


def campaign_terminal_status(root: Path) -> dict[tuple[str, float], str]:
    path = root / "campaign_status_v2.csv"
    out: dict[tuple[str, float], str] = {}
    try:
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                out[(row.get("class", ""), float(row.get("target_psi_deg", "nan")))] = row.get("status", "")
    except (OSError, ValueError):
        pass
    return out


def render(statuses: list[CaseStatus], root: Path, clear: bool, show_active_logs: int) -> None:
    if clear and sys.stdout.isatty():
        print("\033[2J\033[H", end="")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"Mixed-mode FEM/CZM v2.1 monitor — {now}")
    print(f"OUTROOT: {root}")
    print()
    header = (
        f"{'class':<10} {'target':>7} {'state':<10} {'iter':>4} {'alpha':>8} "
        f"{'psi':>8} {'err':>7} {'fit':>4} {'mode':<9} {'step':>11} {'KJ':>8} {'sig':>7} {'age':>5}"
    )
    print(header)
    print("-" * len(header))
    for s in sorted(statuses, key=lambda x: (x.klass.lower(), x.target)):
        step_text = "-"
        if s.step is not None:
            step_text = f"{s.step}/{s.max_steps}" if s.max_steps else str(s.step)
        print(
            f"{s.klass:<10} {s.target:7.1f} {s.state:<10} "
            f"{('-' if s.iteration is None else s.iteration):>4} "
            f"{fmt_num(s.alpha,8,2)} {fmt_num(s.achieved,8,2)} {fmt_num(s.error,7,2)} "
            f"{fmt_bool(s.reliable):>4} {s.mode[:9]:<9} {step_text:>11} "
            f"{fmt_num(s.KJ,8,2)} {fmt_num(s.sig_tip,7,2)} {age_text(s.age_s):>5}"
        )
    counts: dict[str, int] = {}
    for s in statuses:
        counts[s.state] = counts.get(s.state, 0) + 1
    print()
    print("Status: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    completed = sum(s.state in {"CONVERGED", "NOT_CONV"} for s in statuses)
    print(f"Finalized: {completed}/{len(statuses)}")
    if show_active_logs > 0:
        active = [s for s in statuses if s.state in {"RUNNING", "STALE?", "TRIAL_DONE"}]
        if active:
            print("\nRecent active log output:")
            for s in active[:3]:
                log = Path(s.trial_dir) / "run.log"
                lines = read_last_lines(log)
                print(f"\n--- {s.klass} target={s.target:g} iter={s.iteration} ---")
                for line in lines[-show_active_logs:]:
                    print(line)
    print("\nCtrl-C stops only this monitor; it does not stop the campaign.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("runs/mixed_mode_fem_czm_v2_event_controlled_500K"))
    ap.add_argument("--classes", default="ceramic DBTT")
    ap.add_argument("--target-psi-deg", default="-60 -45 -30 -15 0 15 30 45 60")
    ap.add_argument("--interval", type=float, default=15.0)
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--no-clear", action="store_true")
    ap.add_argument("--max-steps-default", type=int, default=3000)
    ap.add_argument("--stale-minutes", type=float, default=15.0)
    ap.add_argument("--show-active-log-lines", type=int, default=0)
    ap.add_argument("--exit-when-final", action="store_true")
    args = ap.parse_args()

    classes = parse_items(args.classes)
    targets = parse_items(args.target_psi_deg, float)
    args.root = args.root.expanduser().resolve()
    expected = [(klass, target) for klass in classes for target in targets]
    if not expected:
        ap.error("No classes or target phase angles were supplied.")

    try:
        while True:
            statuses = [
                status_for_case(args.root, klass, target, args.max_steps_default, args.stale_minutes * 60.0)
                for klass, target in expected
            ]
            terminal = campaign_terminal_status(args.root)
            for s in statuses:
                campaign_state = terminal.get((s.klass, s.target))
                if campaign_state == "failed" and s.state not in {"CONVERGED", "NOT_CONV"}:
                    s.state = "FAILED"
            render(statuses, args.root, not args.no_clear, args.show_active_log_lines)
            all_final = all(s.state in {"CONVERGED", "NOT_CONV", "FAILED"} for s in statuses)
            if args.once or (args.exit_when_final and all_final):
                return 0 if all(s.state != "FAILED" for s in statuses) else 2
            time.sleep(max(1.0, args.interval))
    except KeyboardInterrupt:
        print("\nMonitor stopped. Campaign was not modified.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
