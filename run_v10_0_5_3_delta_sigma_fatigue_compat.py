#!/usr/bin/env python3
"""Audited project-tree launcher for the v10.0.5.3 Delta-sigma campaign.

The launcher performs a construction preflight against the exact current
run_2d/v10.0.3/v10.0.2 source-transform chain before creating a campaign. It
also keeps the campaign CLI synchronized with the inherited parser, streams
child output to both the terminal and its run log, emits periodic heartbeats,
and prints the child log tail on failure.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shlex
import subprocess
import threading
import time

from arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue_audited import (
    validate_source_transform_v10053,
)


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "run_v10_0_5_3_delta_sigma_fatigue.py"
OLD_ENTRY = "arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue"
AUDITED_ENTRY = (
    "arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue_audited"
)


def _load_campaign_module():
    spec = importlib.util.spec_from_file_location("_v10053_delta_sigma_campaign", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load campaign module from {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _remove_option_pair(command: list[str], option: str) -> None:
    while option in command:
        index = command.index(option)
        del command[index:index + 2]


def _replace_entry_module(command: list[str]) -> None:
    try:
        index = command.index(OLD_ENTRY)
    except ValueError as exc:
        raise RuntimeError(
            f"campaign command does not contain expected entry module {OLD_ENTRY}"
        ) from exc
    command[index] = AUDITED_ENTRY


def _stage_from_log(log_path: Path) -> str:
    text = str(log_path)
    if "/calibration/" in text:
        return "calibration"
    if "DeltaSigma_" in text:
        return "fatigue-case"
    return "child-run"


def main() -> int:
    preflight = validate_source_transform_v10053()
    if not preflight.get("source_transform_preflight_passed", False):
        raise RuntimeError("v10.0.5.3 source-transform preflight did not pass")
    print("V10_0_5_3_FATIGUE_SOURCE_PREFLIGHT PASS", flush=True)

    campaign = _load_campaign_module()
    original_base_command = campaign._base_command

    def compatible_base_command(args, outdir, temperature, dU_m):
        command = list(original_base_command(args, outdir, temperature, dU_m))
        _replace_entry_module(command)

        # The active parser has --max-da-per-block-um but no separate target-da
        # option. Cycle-jump resolution remains controlled by target-dB and the
        # finite target-dN limits.
        _remove_option_pair(command, "--target-da-per-block-um")

        # The inherited direct Mode-I v9.11 integration requires competing
        # crystallographic directions whenever anisotropic elasticity is active.
        if "--crystal-compete" not in command:
            try:
                index = command.index("--crystal-aniso") + 1
            except ValueError:
                index = len(command)
            command.insert(index, "--crystal-compete")
        return command

    def verbose_run(command: list[str], env: dict[str, str], log_path: Path) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        child_env = dict(env)
        child_env["PYTHONUNBUFFERED"] = "1"
        stage = _stage_from_log(log_path)
        started = time.monotonic()

        print(f"\nV10_0_5_3 CHILD START stage={stage}", flush=True)
        print(f"log={log_path}", flush=True)
        print(f"command={shlex.join(command)}", flush=True)

        stop_heartbeat = threading.Event()
        with log_path.open("w") as log:
            process = subprocess.Popen(
                command,
                env=child_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            def heartbeat() -> None:
                while not stop_heartbeat.wait(30.0):
                    try:
                        size = log_path.stat().st_size
                    except OSError:
                        size = -1
                    elapsed = time.monotonic() - started
                    print(
                        "V10_0_5_3 HEARTBEAT "
                        f"stage={stage} pid={process.pid} elapsed_s={elapsed:.0f} "
                        f"log_bytes={size}",
                        flush=True,
                    )

            heartbeat_thread = threading.Thread(target=heartbeat, daemon=True)
            heartbeat_thread.start()
            try:
                if process.stdout is None:
                    raise RuntimeError("child process stdout pipe was not created")
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    print(line, end="", flush=True)
                returncode = process.wait()
            finally:
                stop_heartbeat.set()
                heartbeat_thread.join(timeout=1.0)

        elapsed = time.monotonic() - started
        print(
            f"V10_0_5_3 CHILD END stage={stage} pid={process.pid} "
            f"exit={returncode} elapsed_s={elapsed:.1f}",
            flush=True,
        )
        if returncode != 0:
            try:
                lines = log_path.read_text(errors="replace").splitlines()
                tail = "\n".join(lines[-120:])
            except Exception as exc:
                tail = f"<could not read log tail: {type(exc).__name__}: {exc}>"
            raise RuntimeError(
                "command failed with exit code "
                f"{returncode}; see {log_path}\n"
                "----- child log tail -----\n"
                f"{tail}\n"
                "----- end child log tail -----"
            )

    campaign._base_command = compatible_base_command
    campaign._run = verbose_run
    return int(campaign.main() or 0)


if __name__ == "__main__":
    raise SystemExit(main())
