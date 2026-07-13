#!/usr/bin/env python3
"""Shared activation-entropy family for matched Figure 1 Panels C and D.

The same entropy function is used in both calculations:

    S*(sigma,T)/kB = clip[-A_T F_T(T) - A_sigma F_sigma(sigma), S_min, 0]

where

    F_T(T) = h_T(T)/h_T(T0),
    h_T(T) = (T/T_S)^p_T / [1 + (T/T_S)^p_T],

and

    F_sigma(sigma) = x^p_sigma / (1 + x^p_sigma),
    x = |sigma|/sigma_S.

Thus A_T is the magnitude of the temperature-dependent entropy contribution at
T0, while A_sigma is the amplitude of the stress-dependent entropy contribution.
A_sigma is the common waterfall axis in Panels C and D; A_T is the common color
coordinate.

The parent barrier is

    G*(sigma,T) = H(sigma) - T S*(sigma,T) - sigma v(sigma),

using the same rational H(sigma) and v(sigma) forms as the sharp-front model.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numpy as np

KB_J_PER_K = 1.380649e-23
EV_TO_J = 1.602176634e-19


@dataclass(frozen=True)
class EntropyFamily:
    """Two-coordinate entropy family shared by Panels C and D."""

    A_T_kB: float = 4.0
    A_sigma_kB: float = 4.0
    T0_K: float = 300.0
    T_S_K: float = 400.0
    T_gate_power: float = 4.0
    sigma_S_GPa: float = 3.0
    sigma_gate_power: float = 1.0
    S_min_kB: float = -40.0
    S_max_kB: float = 0.0

    @staticmethod
    def _hill(x, power: float):
        x = np.maximum(np.asarray(x, dtype=float), 0.0)
        p = max(float(power), 1e-12)
        xp = np.power(x, p)
        return xp / (1.0 + xp)

    def temperature_gate(self, T_K: float):
        """Saturating temperature gate normalized to one at T0."""
        T = max(float(T_K), 1e-12)
        Ts = max(float(self.T_S_K), 1e-12)
        g = float(self._hill(T / Ts, self.T_gate_power))
        g0 = float(self._hill(max(float(self.T0_K), 1e-12) / Ts, self.T_gate_power))
        return g / max(g0, 1e-30)

    def stress_gate(self, sigma_Pa):
        sigma = np.abs(np.asarray(sigma_Pa, dtype=float))
        s0 = max(float(self.sigma_S_GPa) * 1e9, 1.0)
        return self._hill(sigma / s0, self.sigma_gate_power)

    # Alias retained for simple diagnostics/legacy plotting helpers.
    def gate(self, sigma_Pa):
        return self.stress_gate(sigma_Pa)

    def S_kB(self, sigma_Pa, T_K: float):
        raw = (
            -float(self.A_T_kB) * float(self.temperature_gate(T_K))
            -float(self.A_sigma_kB) * self.stress_gate(sigma_Pa)
        )
        return np.clip(raw, float(self.S_min_kB), float(self.S_max_kB))

    def S_J_per_K(self, sigma_Pa, T_K: float):
        return self.S_kB(sigma_Pa, T_K) * KB_J_PER_K


@dataclass(frozen=True)
class ParentArrheniusBarrier:
    """Common H - T*S - sigma*v barrier used by both Panels C and D."""

    H0_eV: float = 0.8
    sigma0_H_GPa: float = 2.5
    v0_b3: float = 0.6
    sigma0_v_GPa: float = 2.5
    b_m: float = 2.74e-10
    entropy: EntropyFamily = EntropyFamily()
    monotone_stress: bool = True
    monotone_search_GPa: float = 100.0
    monotone_grid_n: int = 4001

    def H_J(self, sigma_Pa):
        sigma = np.abs(np.asarray(sigma_Pa, dtype=float))
        s0 = max(float(self.sigma0_H_GPa) * 1e9, 1.0)
        return float(self.H0_eV) * EV_TO_J / (1.0 + (sigma / s0) ** 2)

    def v_m3(self, sigma_Pa):
        sigma = np.abs(np.asarray(sigma_Pa, dtype=float))
        s0 = max(float(self.sigma0_v_GPa) * 1e9, 1.0)
        v0 = float(self.v0_b3) * float(self.b_m) ** 3
        return v0 / (1.0 + (sigma / s0) ** 1.5)

    def raw_G_J(self, sigma_Pa, T_K: float):
        sigma = np.abs(np.asarray(sigma_Pa, dtype=float))
        G = (
            self.H_J(sigma)
            -float(T_K) * self.entropy.S_J_per_K(sigma, T_K)
            -sigma * self.v_m3(sigma)
        )
        return np.maximum(G, 0.0)

    @lru_cache(maxsize=2048)
    def _sigma_argmin_scalar(self, T_key: float) -> float:
        sg = np.linspace(
            0.0,
            max(float(self.monotone_search_GPa), 1e-6) * 1e9,
            max(int(self.monotone_grid_n), 101),
        )
        G = self.raw_G_J(sg, float(T_key))
        return float(sg[int(np.argmin(G))])

    def G_J(self, sigma_Pa, T_K: float):
        sigma = np.abs(np.asarray(sigma_Pa, dtype=float))
        if self.monotone_stress:
            sstar = self._sigma_argmin_scalar(round(float(T_K), 6))
            sigma = np.minimum(sigma, sstar)
        return self.raw_G_J(sigma, T_K)

    def G_eV(self, sigma_Pa, T_K: float):
        return self.G_J(sigma_Pa, T_K) / EV_TO_J

    def rate_s(self, sigma_Pa, T_K: float, nu0_s: float):
        G = self.G_J(sigma_Pa, T_K)
        expo = -G / max(KB_J_PER_K * float(T_K), 1e-300)
        return float(nu0_s) * np.exp(np.clip(expo, -700.0, 0.0))


# Backward-compatible alias for prior internal tests.
ClassicArrheniusBarrier = ParentArrheniusBarrier


def build_parent_barrier(
    A_T_kB: float,
    A_sigma_kB: float,
    *,
    T0_K: float = 300.0,
    T_S_K: float = 400.0,
    T_gate_power: float = 4.0,
    sigma_S_GPa: float = 3.0,
    sigma_gate_power: float = 1.0,
    S_min_kB: float = -40.0,
    H0_eV: float = 0.8,
    sigma0_H_GPa: float = 2.5,
    v0_b3: float = 0.6,
    sigma0_v_GPa: float = 2.5,
    b_m: float = 2.74e-10,
) -> ParentArrheniusBarrier:
    ent = EntropyFamily(
        A_T_kB=float(A_T_kB),
        A_sigma_kB=float(A_sigma_kB),
        T0_K=float(T0_K),
        T_S_K=float(T_S_K),
        T_gate_power=float(T_gate_power),
        sigma_S_GPa=float(sigma_S_GPa),
        sigma_gate_power=float(sigma_gate_power),
        S_min_kB=float(S_min_kB),
        S_max_kB=0.0,
    )
    return ParentArrheniusBarrier(
        H0_eV=float(H0_eV),
        sigma0_H_GPa=float(sigma0_H_GPa),
        v0_b3=float(v0_b3),
        sigma0_v_GPa=float(sigma0_v_GPa),
        b_m=float(b_m),
        entropy=ent,
        monotone_stress=True,
    )


def multihit_rate(raw_rate_s, m_hits: float = 3.0, tau_s: float = 1e-6):
    """Cooperative renewal rate used for the optional Panel-C crack event."""
    raw = np.asarray(raw_rate_s, dtype=float)
    if float(m_hits) <= 1.0 + 1e-12:
        return raw
    from scipy.special import gammainc
    tau = max(float(tau_s), 1e-30)
    return gammainc(float(m_hits), np.minimum(raw * tau, 1e12)) / tau
