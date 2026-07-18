#!/usr/bin/env python3
"""Streamed launcher for the v10.0.5.5 stochastic VHCF campaign."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shlex
import subprocess
import threading
import time

from arrhenius_fracture.mode_i_first_passage_v10_0_5_5_stochastic_vhcf import (
    validate_source_transform_v10055,
)

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "run_v10_0_5_5_stochastic_vhcf_delta_sigma.py"


def _load_campaign_module():
    spec = importlib.util.spec_from_file_location("_v10055_campaign", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load campaign module from {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    preflight = validate_source_transform_v10055()
    if not preflight.get("v10055_source_transform_preflight_passed", False):
        raise RuntimeError("v10.0.5.5 source-transform preflight failed")
    print("V10_0_5_5_STOCHASTIC_VHCF_SOURCE_PREFLIGHT PASS", flush=True)
    campaign = _load_campaign_module()

    def verbose_run(command: list[str], env: dict[str, str], log_path: Path) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        child_env = dict(env)
        child_env["PYTHONUNBUFFERED"] = "1"
        started = time.monotonic()
        stage = "calibration" if "/calibration/" in str(log_path) else "stochastic-vhcf-case"
        print(f"\nV10_0_5_5 CHILD START stage={stage}", flush=True)
        print(f"log={log_path}", flush=True)
        print(f"command={shlex.join(command)}", flush=True)
        stop = threading.Event()
        with log_path.open("w") as log:
            process = subprocess.Popen(
                command, env=child_env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1,
            )

            def heartbeat():
                while not stop.wait(30.0):
                    size = log_path.stat().st_size if log_path.exists() else -1
                    print(
                        f"V10_0_5_5 HEARTBEAT stage={stage} pid={process.pid} "
                        f"elapsed_s={time.monotonic()-started:.0f} log_bytes={size}",
                        flush=True,
                    )

            thread = threading.Thread(target=heartbeat, daemon=True)
            thread.start()
            try:
                if process.stdout is None:
                    raise RuntimeError("child stdout pipe was not created")
                for line in process.stdout:
                    log.write(line); log.flush(); print(line, end="", flush=True)
                returncode = process.wait()
            finally:
                stop.set(); thread.join(timeout=1.0)
        print(
            f"V10_0_5_5 CHILD END stage={stage} pid={process.pid} "
            f"exit={returncode} elapsed_s={time.monotonic()-started:.1f}",
            flush=True,
        )
        if returncode != 0:
            tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-120:])
            raise RuntimeError(
                f"command failed with exit code {returncode}; see {log_path}\n"
                f"----- child log tail -----\n{tail}\n----- end child log tail -----"
            )

    campaign._base._run = verbose_run
    return int(campaign.main() or 0)


if __name__ == "__main__":
    raise SystemExit(main())
