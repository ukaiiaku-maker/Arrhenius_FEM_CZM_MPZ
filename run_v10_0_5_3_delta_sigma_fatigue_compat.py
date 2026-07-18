#!/usr/bin/env python3
"""macOS/project-tree compatibility launcher for the v10.0.5.3 Δσ campaign.

This launcher changes only campaign command construction and failure reporting:
- removes the retired/undefined ``--target-da-per-block-um`` CLI option;
- supplies ``--crystal-compete``, required by the inherited v9.11 Mode-I path;
- prints the tail of a failed child-run log in the raised exception.

No constitutive, FEM, CZM, MPZ, or fatigue kinetics are changed.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "run_v10_0_5_3_delta_sigma_fatigue.py"


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


def main() -> int:
    campaign = _load_campaign_module()
    original_base_command = campaign._base_command

    def compatible_base_command(args, outdir, temperature, dU_m):
        command = list(original_base_command(args, outdir, temperature, dU_m))

        # The active sharp_front parser has --max-da-per-block-um, but no
        # --target-da-per-block-um. Cycle-jump control remains governed by the
        # inherited target-dB and target-dN limits.
        _remove_option_pair(command, "--target-da-per-block-um")

        # The inherited direct Mode-I v9.11 integration requires both anisotropic
        # elasticity and competing crystallographic directions, even with one
        # active crack front and branching disabled.
        if "--crystal-compete" not in command:
            try:
                index = command.index("--crystal-aniso") + 1
            except ValueError:
                index = len(command)
            command.insert(index, "--crystal-compete")
        return command

    def verbose_run(command: list[str], env: dict[str, str], log_path: Path) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w") as log:
            process = subprocess.run(
                command,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        if process.returncode != 0:
            try:
                lines = log_path.read_text(errors="replace").splitlines()
                tail = "\n".join(lines[-80:])
            except Exception as exc:  # diagnostic path must not mask the run error
                tail = f"<could not read log tail: {type(exc).__name__}: {exc}>"
            raise RuntimeError(
                "command failed with exit code "
                f"{process.returncode}; see {log_path}\n"
                "----- child log tail -----\n"
                f"{tail}\n"
                "----- end child log tail -----"
            )

    campaign._base_command = compatible_base_command
    campaign._run = verbose_run
    return int(campaign.main() or 0)


if __name__ == "__main__":
    raise SystemExit(main())
