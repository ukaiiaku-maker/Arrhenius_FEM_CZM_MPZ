"""FEM/CZM front-engine adapter for the PF v10.2.22 persistent-site MPZ."""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from .config import EV_TO_J, KB
from .mpz_front_engine import MovingProcessZoneFrontEngine
from .persistent_site_registry_v100514 import PersistentSiteRowV100514
from .persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
)

MODEL_ID = "FEM_CZM_persistent_site_front_engine_v10_0_5_14"


class PersistentSiteMovingProcessZoneFrontEngineV100514(MovingProcessZoneFrontEngine):
    """Production front engine with trial/commit persistent signed MPZ state."""

    persistent_site_source_active = True
    state_model = "moving_pz"
    _candidate_default: PersistentSiteRowV100514 | None = None
    _kernel_default: SignedShieldingKernelV100514 | None = None

    @classmethod
    def configure(
        cls,
        candidate: PersistentSiteRowV100514,
        kernel: SignedShieldingKernelV100514,
    ) -> None:
        cls._candidate_default = copy.deepcopy(candidate.validate())
        cls._kernel_default = copy.deepcopy(kernel)

    @classmethod
    def clear_configuration(cls) -> None:
        cls._candidate_default = None
        cls._kernel_default = None

    def reset(self):
        candidate = type(self)._candidate_default
        kernel = type(self)._kernel_default
        if candidate is None or kernel is None:
            raise RuntimeError("v10.0.5.14 persistent-site engine was not configured")
        blunting_length = float(getattr(self.mpz_config, "blunting_length_m", 0.5e-6))
        max_cfl = float(getattr(self.mpz_config, "max_transport_cfl", 0.35))
        max_substeps = int(getattr(self.mpz_config, "max_transport_substeps", 2000))
        self.mpz_state = PersistentSiteSignedMPZStateV100514(
            candidate,
            kernel,
            G_Pa=self.G,
            nu=self.nu,
            b_m=self.b,
            r0_m=self.f.r0,
            blunting_length_m=blunting_length,
            wake_shielding=False,
            max_transport_cfl=max_cfl,
            max_transport_substeps=max_substeps,
        )
        self.N_em = 0.0
        self.B = 0.0
        self.a_adv = 0.0
        self.n_adv = 0
        self.W_emit = 0.0
        self.t = 0.0
        self.K_prev = None
        self._lambda_c_prev = None
        self._K_cleave_prev = None
        self._last_pre_renewal_state = None
        self.accepted_constitutive_steps = 0
        self.rejected_constitutive_trials = 0

    def sigma_back(self) -> float:
        _, sigma = self.mpz_state.backstress()
        return float(np.mean(sigma))

    def _two_channel_drive(self) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
        context = getattr(self, "_mm", None)
        latest = {} if context is None else dict(getattr(context, "latest", {}) or {})
        reliable = bool(latest.get("two_channel_drive_reliable", False))
        factors = np.asarray(
            latest.get("two_channel_drive_factors", []), dtype=float
        ).reshape(-1)
        tau = np.asarray(
            latest.get("two_channel_tau_signed_Pa", []), dtype=float
        ).reshape(-1)
        if not reliable or factors.shape != (2,) or tau.shape != (2,):
            raise RuntimeError(
                "v10.0.5.14 persistent signed emission requires a reliable "
                "two-channel FEM tensor drive"
            )
        return factors, tau, latest

    def _emission_rate_per_site(self, stress_Pa: float, T_K: float) -> float:
        stress = max(float(stress_Pa), 0.0)
        barrier = float(self.eb.G_barrier(np.asarray([stress]), T_K, self.b)[0])
        return float(
            max(float(self.f.nu0_e), 0.0)
            * math.exp(
                float(
                    np.clip(
                        -barrier / max(KB * float(T_K), 1.0e-30),
                        -700.0,
                        0.0,
                    )
                )
            )
        )

    def _trial_constitutive_update(
        self,
        *,
        dt_s: float,
        T_K: float,
        opening_stress_Pa: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> tuple[PersistentSiteSignedMPZStateV100514, dict[str, Any]]:
        trial = self.mpz_state.copy()
        kinetics = trial.evolve(
            dt_s=dt_s,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
            emission_rate_function=self._emission_rate_per_site,
        )
        return trial, kinetics

    def step_drives(
        self,
        K_cleave,
        K_emit,
        T,
        dt,
        metadata: dict[str, Any] | None = None,
    ):
        dt = max(float(dt), 0.0)
        factors, tau_signed, drive_metadata = self._two_channel_drive()
        Ksh0 = self.K_shield()
        r0 = self.r_eff()
        opening_stress = max(float(K_emit) - Ksh0, 0.0) / math.sqrt(
            2.0 * math.pi * max(r0, 1.0e-30)
        )
        sigma_cap_active = bool(
            self.f.sigma_cap > 0.0 and opening_stress > self.f.sigma_cap
        )
        if self.f.sigma_cap > 0.0:
            opening_stress = min(opening_stress, self.f.sigma_cap)

        # Trial/commit: no emission survives a failed constitutive evaluation.
        try:
            trial, kinetics = self._trial_constitutive_update(
                dt_s=dt,
                T_K=float(T),
                opening_stress_Pa=opening_stress,
                drive_factors=factors,
                tau_signed_Pa=tau_signed,
            )
        except BaseException:
            self.rejected_constitutive_trials += 1
            raise
        self.mpz_state = trial
        self.accepted_constitutive_steps += 1
        self.W_emit += (
            opening_stress
            * self.b
            * self.f.L_pz
            * float(kinetics.get("dN_emit", 0.0))
        )
        self._sync_compat()

        sig_cleave = self.sigma_tip(K_cleave)
        lam_c, lam_c_raw, Gc = self.lambda_cleave(sig_cleave, T)
        if self.f.tau_B > 0.0 and dt > 0.0:
            self.B *= math.exp(-min(dt / self.f.tau_B, 80.0))
        leff = self._logmean_rate(self._lambda_c_prev, lam_c)
        self.B += leff * dt
        self._lambda_c_prev = lam_c
        self._K_cleave_prev = float(K_cleave)
        self.K_prev = float(K_cleave)
        self.t += dt

        renew = self._renew(dt)
        state_diag = self.mpz_state.diagnostics(
            self.G, self.nu, self.b, self.f.r0, self.f.c_blunt
        )
        last_emission = self.mpz_state.last_emission
        aggregate = np.asarray(
            last_emission.get("aggregate_hazard_initial_by_system_s", np.zeros(2)),
            dtype=float,
        )
        out = {
            **renew,
            "sigma_tip": float(sig_cleave),
            "sigma_emit_tip": float(opening_stress),
            "sigma_back": float(self.sigma_back()),
            "lambda_e": float(np.sum(aggregate)),
            "lambda_c": float(lam_c),
            "lambda_c_raw": float(lam_c_raw),
            "B": float(self.B),
            "N_em": float(self.N_em),
            "r_eff": float(self.r_eff()),
            "dG_emb_eV": 0.0,
            "G_cleave_eff_eV": float(Gc / EV_TO_J),
            **self.cleavage_diagnostics(sig_cleave, T),
            "G_emit_eV": float(
                self.eb.G_barrier(
                    np.asarray([max(opening_stress, 0.0)]), T, self.b
                )[0]
                / EV_TO_J
            ),
            "W_emit": float(self.W_emit),
            "sigma_tip_uncapped": float(
                max(float(K_emit) - Ksh0, 0.0)
                / math.sqrt(2.0 * math.pi * max(r0, 1.0e-30))
            ),
            "sigma_cap_active": sigma_cap_active,
            "dN_emit_raw": float(kinetics.get("dN_emit", 0.0)),
            "dN_cap_active": False,
            "N_sat_factor": 1.0,
            "N_sat_active": False,
            "front_state_model_code": 2.0,
            "persistent_site_source_active": True,
            "persistent_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "available_site_fraction": 1.0,
            "source_sites_refreshed": 0.0,
            "trial_commit_state_active": True,
            "accepted_constitutive_steps": self.accepted_constitutive_steps,
            "rejected_constitutive_trials": self.rejected_constitutive_trials,
            "two_channel_drive_factors": factors.tolist(),
            "two_channel_tau_signed_Pa": tau_signed.tolist(),
            "two_channel_drive_reliable": True,
            "two_channel_names": drive_metadata.get("two_channel_names"),
            **kinetics,
            **state_diag,
        }
        if metadata:
            out.update(metadata)
        return out

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        candidate = cls._candidate_default
        kernel = cls._kernel_default
        return {
            "model_id": MODEL_ID,
            "persistent_site_source": True,
            "finite_source_inventory": False,
            "source_depletion_on_emission": False,
            "source_refresh_on_crack_advance": False,
            "front_width_grid_independent": True,
            "emission_integrator": "implicit_backward_euler_backstress_complementarity",
            "backstress_population": "unsigned_mobile_plus_retained",
            "shielding_population": "signed_retained_only",
            "trial_commit_state": True,
            "candidate": None if candidate is None else candidate.candidate_id,
            "candidate_fingerprint": (
                None if candidate is None else candidate.fingerprint()
            ),
            "kernel_source": None if kernel is None else kernel.source_path,
        }


__all__ = [
    "MODEL_ID",
    "PersistentSiteMovingProcessZoneFrontEngineV100514",
]
