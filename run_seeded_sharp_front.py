#!/usr/bin/env python3
"""Run the active sharp-front/CZM solver with a seeded Poisson first-passage clock.

All constitutive parameters remain unchanged. ``--solver-seed`` initializes the
random-number stream used to realize the nonhomogeneous Poisson cleavage process.
For each renewal, the cumulative physical hazard action is compared with

    H_n = -log(U_n),  U_n ~ Uniform(0,1),

so H_n ~ Exp(1). These thresholds are latent random variates implied by the
hazard model, not user-selected material parameters.

The wrapper avoids modifying installed package files. Internally it normalizes
the existing deterministic B>=1 clock by the current Exp(1) variate. The physical
cleavage rate reported in diagnostics is restored before returning from step().
"""
from __future__ import annotations

import argparse
import inspect
import json
from pathlib import Path
import sys

import numpy as np


def _out_dir_from_args(argv: list[str]) -> Path | None:
    for i, token in enumerate(argv[:-1]):
        if token == "--out":
            return Path(argv[i + 1])
    return None


def install_seeded_first_passage(seed: int) -> None:
    import arrhenius_fracture.sharp_front as sf

    source = inspect.getsource(sf.FrontEngine.step)
    if "max_advances_per_step" not in source:
        raise RuntimeError(
            "The seeded wrapper requires the adaptive-CZM FrontEngine that "
            "retains excess renewal action and limits topology insertion to one "
            "advance per mechanics solve. Apply the current adaptive-CZM patches."
        )

    # A dedicated stream isolates physical first-passage randomness from any
    # unrelated NumPy draws made elsewhere in the solver.
    ss = np.random.SeedSequence([int(seed), 0xC1EA6E])
    rng = np.random.default_rng(ss)
    tiny = float(np.finfo(float).tiny)

    original_reset = sf.FrontEngine.reset
    original_lambda_cleave = sf.FrontEngine.lambda_cleave
    original_step = sf.FrontEngine.step

    def draw_action() -> float:
        # Equivalent to -log(U), but rng.exponential is numerically robust.
        return float(max(rng.exponential(1.0), tiny))

    def reset(self):
        original_reset(self)
        self.first_passage_action_target = draw_action()
        self.first_passage_event_count = 0

    def lambda_cleave(self, sig_tip, T):
        lam_phys, lam_raw, geff = original_lambda_cleave(self, sig_tip, T)
        target = float(max(getattr(self, "first_passage_action_target", 1.0), tiny))
        # Existing B is retained as normalized progress B=H_residual/H_target.
        return lam_phys / target, lam_raw, geff

    def step(self, K, T, dt):
        target_old = float(max(
            getattr(self, "first_passage_action_target", 1.0), tiny
        ))
        info = original_step(self, K, T, dt)

        # Restore the physical rate in output. The base step integrated the
        # normalized rate lambda_phys/H_target solely to reuse its B>=1 logic.
        lambda_normalized = float(info.get("lambda_c", 0.0))
        info["lambda_c_normalized_per_s"] = lambda_normalized
        info["lambda_c"] = lambda_normalized * target_old

        n_fire = int(info.get("n_fire", 0) or 0)
        if n_fire > 1:
            raise RuntimeError(
                "More than one topology event was accepted in one mechanics "
                "solve. Seeded sequential renewal requires the adaptive-CZM "
                "one-event cap."
            )

        if n_fire == 1:
            # Base step subtracted one in old-target-normalized units. Convert
            # remaining action back to physical H, draw the next model-implied
            # Exp(1) target, then renormalize for the next accepted step.
            residual_action = max(float(self.B), 0.0) * target_old
            target_new = draw_action()
            self.first_passage_action_target = target_new
            self.first_passage_event_count = int(
                getattr(self, "first_passage_event_count", 0)
            ) + 1
            self.B = residual_action / target_new
            info["B"] = float(self.B)

        target_now = float(max(
            getattr(self, "first_passage_action_target", target_old), tiny
        ))
        info["first_passage_action_target"] = target_now
        info["first_passage_action_residual"] = float(max(self.B, 0.0) * target_now)
        info["first_passage_event_count"] = int(
            getattr(self, "first_passage_event_count", 0)
        )
        return info

    sf.FrontEngine.reset = reset
    sf.FrontEngine.lambda_cleave = lambda_cleave
    sf.FrontEngine.step = step


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--solver-seed", type=int, required=True)
    known, remaining = parser.parse_known_args()

    install_seeded_first_passage(known.solver_seed)

    out = _out_dir_from_args(remaining)
    if out is not None:
        out.mkdir(parents=True, exist_ok=True)
        (out / "seeded_first_passage_config.json").write_text(
            json.dumps(
                {
                    "solver_seed": int(known.solver_seed),
                    "random_process": "nonhomogeneous Poisson cleavage renewal",
                    "latent_action_distribution": "H_n=-ln(U_n), U_n~Uniform(0,1); H_n~Exp(1)",
                    "constitutive_parameters_changed_between_replicates": False,
                    "mesh_seed_changed_between_replicates": False,
                    "purpose": "physical Arrhenius first-passage scatter",
                },
                indent=2,
            )
        )

    import arrhenius_fracture.sharp_front as sf
    sys.argv = [sys.argv[0], *remaining]
    sf.main()


if __name__ == "__main__":
    main()
