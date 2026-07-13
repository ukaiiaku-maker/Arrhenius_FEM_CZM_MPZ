"""Fully Arrhenius emission--Peierls--Taylor plastic-event chain for S-N models.

The S-N initiation solvers use the same scaled EXP-floor free-energy families as
``fatigue_v1.py``.  No athermal Taylor stress, Peierls floor, or quasi-static
return surface is inserted.  The elementary event rates are evaluated directly:

    lambda_e = nu_e exp[-DG_e(sigma,T)/(kBT)]
    lambda_P = nu_P exp[-DG_P(sigma,T)/(kBT)]
    lambda_T = nu_T exp[-DG_T(phi_T(rho)*sigma,T)/(kBT)]

with the Taylor node amplification

    phi_T(rho) = min[1/(2 b sqrt(rho)), phi_max].

Peierls and Taylor are sequential glide/depinning resistances and therefore
combine by reciprocal residence times.  Emission is placed in series with that
mobility branch for a completed plastic-flow event.  This is a rate-space
construction; stresses are never added and no hard yield gate is imposed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict
import numpy as np

from .fatigue_v1 import ExpFloorBarrierParams, ScaledExpFloorBarrier


@dataclass
class ArrheniusPlasticChain:
    emit: ScaledExpFloorBarrier
    peierls: ScaledExpFloorBarrier
    taylor: ScaledExpFloorBarrier
    b_m: float
    phi_taylor_max: float = 20.0
    plastic_event_strain: float = 1.0e-5

    @staticmethod
    def _series_rate(*rates):
        inv = None
        for r in rates:
            a = np.maximum(np.asarray(r, dtype=float), 1.0e-300)
            inv = 1.0/a if inv is None else inv + 1.0/a
        return 1.0 / np.maximum(inv, 1.0e-300)

    def taylor_phi(self, rho_m2):
        rho = np.maximum(np.asarray(rho_m2, dtype=float), 1.0e6)
        delta = 1.0 / (2.0*np.sqrt(rho))
        return np.minimum(delta / max(self.b_m, 1e-30), max(float(self.phi_taylor_max), 1.0))

    def rates(self, sigma_eq_Pa, rho_m2, T_K: float) -> Dict[str, np.ndarray]:
        sig = np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
        phi_T = self.taylor_phi(rho_m2)
        lam_e = self.emit.rate(sig, T_K)
        lam_P = self.peierls.rate(sig, T_K)
        lam_T = self.taylor.rate(phi_T * sig, T_K)
        lam_escape = self._series_rate(lam_P, lam_T)
        lam_flow = self._series_rate(lam_e, lam_escape)
        dot_ep = max(float(self.plastic_event_strain), 0.0) * lam_flow
        return {
            "lambda_emit": lam_e,
            "lambda_peierls": lam_P,
            "lambda_taylor": lam_T,
            "lambda_escape": lam_escape,
            "lambda_flow": lam_flow,
            "dot_ep": dot_ep,
            "phi_taylor": phi_T,
        }

    def cycle_integrals(self, sigma_hist_Pa, rho_m2, T_K: float, frequency_Hz: float):
        """Integrate the Arrhenius chain over one representative cycle.

        The original implementation looped over phase points in Python.  The
        barrier kernels are NumPy-vectorized, so evaluating the complete
        ``(n_phase, n_state)`` array at once is algebraically identical and
        substantially faster for the large barrier-correlation sweeps.
        """
        sig_hist = np.asarray(sigma_hist_Pa, dtype=float)
        if sig_hist.ndim == 1:
            sig_hist = sig_hist[:, None]
        nphase = max(sig_hist.shape[0], 1)
        r = self.rates(sig_hist, rho_m2, T_K)
        period = 1.0 / max(float(frequency_Hz), 1e-300)
        out = {
            "mu_" + k.replace("lambda_", ""): np.mean(np.asarray(v, dtype=float), axis=0) * period
            for k, v in r.items() if k.startswith("lambda_")
        }
        out["dep_eq_per_cycle"] = np.mean(np.asarray(r["dot_ep"], dtype=float), axis=0) * period
        # phi_T is state-dependent but not phase-dependent.
        out["phi_taylor_mean"] = np.asarray(r["phi_taylor"], dtype=float)
        return out

    def barrier_diagnostics(self, sigma_eq_Pa, rho_m2, T_K: float):
        sig = np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
        phi = self.taylor_phi(rho_m2)
        return {
            "G_emit_eV": self.emit.deltaG_eV(sig, T_K),
            "G_peierls_eV": self.peierls.deltaG_eV(sig, T_K),
            "G_taylor_eV": self.taylor.deltaG_eV(phi*sig, T_K),
            "phi_taylor": phi,
        }


def build_chain_from_namespace(args, b_m: float) -> ArrheniusPlasticChain:
    """Build the S-N plastic-event chain from argparse-like attributes.

    Defaults reproduce the selected case-64-M1 fatigue scaling:
      emission energy/entropy scale 0.75/0.75,
      Peierls 0.00375/0.00375,
      Taylor  0.015/0.015.

    The general fatigue-model values remain CLI-overridable.
    """
    base = ExpFloorBarrierParams.preset(getattr(args, "exp_system", "W[100]"))

    # Representative-map and inverse-design studies override the complete
    # EXP-floor free-energy surface.  Earlier S-N map drivers populated these
    # namespace fields, but the chain builder only consumed ``a`` and ``n``;
    # consequently G00, gT, sigc0, sT, Tref, and floor-fraction sweeps were not
    # actually reaching the plastic hazards.  Apply every supported override
    # explicitly and with units made unambiguous here.
    override_map = {
        "exp_G00_eV": ("G00_eV", 1.0),
        "exp_gT_eV_per_K": ("gT_eV_per_K", 1.0),
        "exp_sigc0_Pa": ("sigc0_Pa", 1.0),
        "exp_sigc0_GPa": ("sigc0_Pa", 1.0e9),
        "exp_sT_Pa_per_K": ("sT_Pa_per_K", 1.0),
        "exp_sT_MPa_per_K": ("sT_Pa_per_K", 1.0e6),
        "exp_Tref_K": ("Tref_K", 1.0),
        "exp_a": ("a", 1.0),
        "exp_n": ("n", 1.0),
        "exp_floor_frac": ("Gfloor_fraction", 1.0),
        "exp_Gfloor_fraction": ("Gfloor_fraction", 1.0),
        "exp_Gfloor_min_eV": ("Gfloor_min_eV", 1.0),
        "exp_Gfloor_max_fraction": ("Gfloor_max_fraction", 1.0),
    }
    for arg_name, (field_name, scale) in override_map.items():
        value = getattr(args, arg_name, None)
        if value is not None:
            setattr(base, field_name, float(value) * scale)

    emit = ScaledExpFloorBarrier(
        base=base, mechanism="surface_dislocation_emission",
        energy_scale=float(getattr(args, "emit_energy_scale", 0.75)),
        entropy_scale=float(getattr(args, "emit_entropy_scale", 0.75)),
        stress_scale=float(getattr(args, "emit_stress_scale", 1.0)),
        rate_prefactor=float(getattr(args, "nu0_emit_pz", 1.0e11)),
    )
    peierls = ScaledExpFloorBarrier(
        base=base, mechanism="peierls_glide",
        energy_scale=float(getattr(args, "peierls_energy_scale", 0.00375)),
        entropy_scale=float(getattr(args, "peierls_entropy_scale", 0.00375)),
        stress_scale=float(getattr(args, "peierls_stress_scale", 1.0)),
        rate_prefactor=float(getattr(args, "nu0_peierls", 1.0e11)),
    )
    taylor = ScaledExpFloorBarrier(
        base=base, mechanism="taylor_junction_depinning",
        energy_scale=float(getattr(args, "taylor_energy_scale", 0.015)),
        entropy_scale=float(getattr(args, "taylor_entropy_scale", 0.015)),
        stress_scale=float(getattr(args, "taylor_stress_scale", 1.0)),
        rate_prefactor=float(getattr(args, "nu0_taylor", 1.0e11)),
    )
    return ArrheniusPlasticChain(
        emit=emit, peierls=peierls, taylor=taylor, b_m=float(b_m),
        phi_taylor_max=float(getattr(args, "phi_taylor_max", 20.0)),
        plastic_event_strain=float(getattr(args, "plastic_event_strain", 1.0e-5)),
    )
