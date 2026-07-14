"""FrontEngine implementation using the moving process-zone state."""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from .config import KB, EV_TO_J
from .moving_process_zone import MovingProcessZoneConfig, MovingProcessZoneState
from .sharp_front import FrontEngine


class MovingProcessZoneFrontEngine(FrontEngine):
    """Drop-in sharp-front engine with finite sites and direct K shielding.

    Geometry, J integration, branching, anisotropy and crack backends remain in
    :mod:`sharp_front`.  This class replaces only the front-local scalar closure.
    """

    state_model = "moving_pz"

    def __init__(self, fcfg, cleave_barrier, emit_barrier, G_shear, nu, b,
                 mpz_config: MovingProcessZoneConfig):
        self.mpz_config = copy.deepcopy(mpz_config)
        self.mpz_state: MovingProcessZoneState | None = None
        self._last_pre_renewal_state: MovingProcessZoneState | None = None
        self._lambda_c_prev: float | None = None
        self._K_cleave_prev: float | None = None
        super().__init__(fcfg, cleave_barrier, emit_barrier, G_shear, nu, b)

    def reset(self):
        self.mpz_state = MovingProcessZoneState(self.mpz_config)
        self.N_em = 0.0  # compatibility projection: retained count only
        self.B = 0.0
        self.a_adv = 0.0
        self.n_adv = 0
        self.W_emit = 0.0
        self.t = 0.0
        self.K_prev = None
        self._lambda_c_prev = None
        self._K_cleave_prev = None
        self._last_pre_renewal_state = None

    def _sync_compat(self) -> None:
        self.N_em = float(self.mpz_state.retained_count)

    def clone_split(self, daughter_fraction=0.5):
        frac = float(np.clip(daughter_fraction, 0.0, 1.0))
        child = copy.deepcopy(self)
        child.mpz_state = self.mpz_state.split(frac)
        B0, W0 = float(self.B), float(self.W_emit)
        child.B = B0 * frac
        child.W_emit = W0 * frac
        child.a_adv = 0.0
        child.n_adv = 0
        child._last_pre_renewal_state = None
        self.B = B0 * (1.0 - frac)
        self.W_emit = W0 * (1.0 - frac)
        self._sync_compat(); child._sync_compat()
        return child

    # ------------------------------------------------------------------
    # Direct process-zone mechanics
    # ------------------------------------------------------------------
    def K_shield(self) -> float:
        return self.mpz_state.shielding_K(self.G, self.nu, self.b)

    def r_eff(self):
        return self.mpz_state.blunted_radius(self.f.r0, self.f.c_blunt, self.b)

    def sigma_tip(self, K):
        K_eff = max(float(K) - self.K_shield(), 0.0)
        s = K_eff / np.sqrt(2.0 * np.pi * max(self.r_eff(), 1.0e-30))
        if self.f.sigma_cap > 0:
            s = min(s, self.f.sigma_cap)
        return float(s)

    def sigma_back(self):
        """Compatibility diagnostic: K_sh converted to an equivalent tip stress.

        This value is never subtracted again in the constitutive equations.
        """
        return float(self.K_shield() / np.sqrt(2.0 * np.pi * max(self.r_eff(), 1.0e-30)))

    def e_stored(self):
        rho = self.f.rho0 + self.mpz_state.retained_count / max(self.f.L_pz ** 2, 1.0e-30)
        return 0.5 * self.G * self.b ** 2 * rho

    def dG_emb(self):
        # Deliberately disabled in v9.  Stored-energy-assisted cleavage may be
        # reintroduced only through an independently constrained physical law.
        return 0.0

    def lambda_emit(self, sig_tip, T):
        s_eff = max(float(sig_tip), 0.0)
        Gstar = float(self.eb.G_barrier(np.array([s_eff]), T, self.b)[0])
        x = -Gstar / max(KB * T, 1.0e-30)
        return self.f.nu0_e * np.exp(np.clip(x, -700.0, 0.0)), s_eff, Gstar

    def lambda_cleave(self, sig_tip, T):
        s_eff = max(float(sig_tip), 0.0)
        Gstar = float(self.cb.G_barrier(np.array([s_eff]), T, self.b)[0])
        x = -Gstar / max(KB * T, 1.0e-30)
        lam_raw = self.f.nu0_c * np.exp(np.clip(x, -700.0, 0.0))
        m = max(self.f.m_hits, 1.0)
        if m > 1.0 + 1.0e-12:
            from scipy.special import gammainc
            tau = max(self.f.tau_c, 1.0e-30)
            lam = gammainc(m, min(lam_raw * tau, 1.0e12)) / tau
        else:
            lam = lam_raw
        return float(lam), float(lam_raw), float(Gstar)

    def cleavage_diagnostics(self, sig_tip, T):
        s_eff = max(float(sig_tip), 0.0)
        d = self.cb.diagnostics(np.array([s_eff]), T, self.b)
        Gstar = float(d["G_eV"][0] * EV_TO_J)
        return {
            "sigma_cleave_eff_Pa": s_eff,
            "G_cleave_raw_eV": Gstar / EV_TO_J,
            "G_cleave_eff_eV": Gstar / EV_TO_J,
            "S_cleave_kB": float(d["S_kB"][0]),
            "dGcleave_dsigma_eV_per_GPa": float(d["dG_dsigma_eV_per_GPa"][0]),
            "vstar_cleave_b3": float(d["vstar_b3"][0]),
            "cleave_barrier_kind_code": 1.0 if str(getattr(self.cb, "barrier_kind", "classic")) == "exp_floor" else 0.0,
        }

    @staticmethod
    def _logmean_rate(lam1: float | None, lam2: float) -> float:
        if lam1 is None:
            return max(float(lam2), 0.0)
        lo, hi = sorted((max(float(lam1), 0.0), max(float(lam2), 0.0)))
        if lo <= 0.0:
            return 0.5 * hi
        if abs(hi - lo) <= 1.0e-12 * max(hi, 1.0e-300):
            return hi
        return (hi - lo) / math.log(hi / lo)

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        sig = self.sigma_tip(K_cleave)
        lam, _, _ = self.lambda_cleave(sig, T)
        return float(max(self._logmean_rate(self._lambda_c_prev, lam) * max(float(dt), 0.0), 0.0))

    def predict_clock_increment(self, K, T, dt):
        return self.predict_clock_increment_drives(K, K, T, dt)

    def _renew(self, dt: float) -> dict[str, Any]:
        self._sync_compat()
        Npre = float(self.N_em)
        Ksh_pre = self.K_shield()
        rpre = self.r_eff()
        mobile_pre = float(self.mpz_state.mobile_count)
        retained_pre = float(self.mpz_state.retained_count)
        site_fraction_pre = float(self.mpz_state.available_site_fraction)
        local_slip_pre = float(self.mpz_state.local_slip_count())
        emitted_total_pre = float(self.mpz_state.emitted_total)
        self._last_pre_renewal_state = self.mpz_state.copy()

        if not np.isfinite(self.B):
            self.B = 0.0
        navail = int(np.floor(min(max(self.B, 0.0), 1.0e7)))
        max_fire = float(getattr(self.f, "max_advances_per_step", float("inf")))
        nfire = min(navail, max(int(max_fire), 0)) if np.isfinite(max_fire) else navail
        fired = nfire >= 1
        wake = {"wake_mobile": 0.0, "wake_retained": 0.0, "wake_slip": 0.0,
                "source_sites_refreshed": 0.0}
        if fired:
            self.B -= nfire
            distance = self.f.da * nfire
            wake = self.mpz_state.advance(distance)
            self.a_adv += distance
            self.n_adv += nfire
        self._sync_compat()
        return {
            "fired": bool(fired), "n_fire": int(nfire),
            "n_fire_available": int(navail),
            "v_crack": self.f.da * nfire / dt if dt > 0 else 0.0,
            "N_em_pre_renewal": Npre,
            "N_em_retained": float(self.N_em),
            "N_em_shed_to_wake": float(wake["wake_retained"]),
            "sigma_back_pre_renewal": float(Ksh_pre / np.sqrt(2.0 * np.pi * max(rpre, 1.0e-30))),
            "r_eff_pre_renewal": float(rpre),
            "mpz_K_shield_pre_renewal_Pa_sqrt_m": float(Ksh_pre),
            "mpz_mobile_pre_renewal": mobile_pre,
            "mpz_retained_pre_renewal": retained_pre,
            "mpz_available_site_fraction_pre_renewal": site_fraction_pre,
            "mpz_local_slip_pre_renewal": local_slip_pre,
            "mpz_emitted_total_pre_renewal": emitted_total_pre,
            "dG_emb_pre_renewal_eV": 0.0,
            "mpz_wake_mobile_block": float(wake["wake_mobile"]),
            "mpz_wake_retained_block": float(wake["wake_retained"]),
            "mpz_wake_slip_block": float(wake["wake_slip"]),
            "mpz_source_sites_refreshed_on_advance": float(wake["source_sites_refreshed"]),
        }

    def restore_geometry_veto(self, n_restore: int) -> None:
        """Undo the last moving-frame renewal when geometry rejected the event."""
        if self._last_pre_renewal_state is not None:
            self.mpz_state = self._last_pre_renewal_state.copy()
        self.B += float(max(int(n_restore), 0))
        self.a_adv = max(self.a_adv - self.f.da * max(int(n_restore), 0), 0.0)
        self.n_adv = max(self.n_adv - max(int(n_restore), 0), 0)
        self._sync_compat()

    def step_drives(self, K_cleave, K_emit, T, dt, metadata: dict[str, Any] | None = None):
        dt = max(float(dt), 0.0)
        sig_emit = self.sigma_tip(K_emit)
        Ksh0 = self.K_shield()
        r0 = self.r_eff()
        sigma_uncapped = max(float(K_emit) - Ksh0, 0.0) / np.sqrt(2.0 * np.pi * max(r0, 1.0e-30))
        sigma_cap_active = bool(self.f.sigma_cap > 0 and sigma_uncapped > self.f.sigma_cap)

        lam_e, sig_em_eff, Ge = self.lambda_emit(sig_emit, T)
        kinetics = self.mpz_state.evolve(
            dt, T, sig_em_eff, self.b,
            emission_hazard_integral=lam_e * dt,
        )
        self.W_emit += sig_em_eff * self.b * self.f.L_pz * kinetics["dN_emit"]
        self._sync_compat()

        sig_cleave = self.sigma_tip(K_cleave)
        lam_c, lam_c_raw, Gc = self.lambda_cleave(sig_cleave, T)
        if self.f.tau_B > 0 and dt > 0:
            self.B *= np.exp(-min(dt / self.f.tau_B, 80.0))
        leff = self._logmean_rate(self._lambda_c_prev, lam_c)
        self.B += leff * dt
        self._lambda_c_prev = lam_c
        self._K_cleave_prev = float(K_cleave)
        self.K_prev = float(K_cleave)
        self.t += dt

        renew = self._renew(dt)
        mpzd = self.mpz_state.diagnostics(self.G, self.nu, self.b, self.f.r0, self.f.c_blunt)
        out = {
            **renew,
            "sigma_tip": float(sig_cleave),
            "sigma_emit_tip": float(sig_emit),
            "sigma_back": float(self.sigma_back()),
            "lambda_e": float(lam_e), "lambda_c": float(lam_c),
            "lambda_c_raw": float(lam_c_raw), "B": float(self.B),
            "N_em": float(self.N_em), "r_eff": float(self.r_eff()),
            "dG_emb_eV": 0.0, "G_cleave_eff_eV": float(Gc / EV_TO_J),
            **self.cleavage_diagnostics(sig_cleave, T),
            "G_emit_eV": float(Ge / EV_TO_J), "W_emit": float(self.W_emit),
            "sigma_tip_uncapped": float(sigma_uncapped),
            "sigma_cap_active": sigma_cap_active,
            "dN_emit_raw": float(lam_e * dt),
            "dN_cap_active": False,
            "N_sat_factor": 1.0, "N_sat_active": False,
            "front_state_model_code": 1.0,
            **kinetics, **mpzd,
        }
        if metadata:
            out.update(metadata)
        return out

    def step(self, K, T, dt):
        return self.step_drives(K, K, T, dt)

    # ------------------------------------------------------------------
    # Unified cyclic loading using the same emission/transport parameters
    # ------------------------------------------------------------------
    def predict_fatigue_cycle(self, waveform, T_K: float, n_phase: int = 96) -> dict[str, float]:
        n = max(int(n_phase), 8)
        phase = (np.arange(n, dtype=float) + 0.5) * 2.0 * np.pi / n
        Kvals = waveform.K_phase(phase)
        dtp = waveform.period_s / n
        sig = np.array([self.sigma_tip(float(k)) for k in Kvals], dtype=float)
        lam_e = np.array([self.lambda_emit(float(s), T_K)[0] for s in sig], dtype=float)
        H_emit = float(np.sum(lam_e) * dtp)
        p_emit = 1.0 - math.exp(-min(max(H_emit, 0.0), 700.0))
        dN_emit = float(np.sum(self.mpz_state.available_sites) * p_emit)
        lam_c = np.array([self.lambda_cleave(float(s), T_K)[0] for s in sig], dtype=float)
        mu_c = float(np.sum(lam_c) * dtp)
        w = np.maximum(lam_e, 0.0)
        avg_sig = float(np.sum(w * sig) / np.sum(w)) if np.sum(w) > 0 else float(np.mean(sig))
        return {
            "H_emit_per_cycle": H_emit,
            "dN_emit_per_cycle": dN_emit,
            "mu_cleave_per_cycle": mu_c,
            "avg_sigma_emit_eff": avg_sig,
            "avg_sigma_tip": float(np.mean(sig)),
            "max_sigma_tip": float(np.max(sig)),
        }

    def commit_fatigue_block(self, waveform, T_K: float, cycles: float,
                             n_phase: int = 96) -> dict[str, Any]:
        cycles = max(float(cycles), 0.0)
        pred = self.predict_fatigue_cycle(waveform, T_K, n_phase)
        dt_block = cycles * waveform.period_s
        kinetics = self.mpz_state.evolve(
            dt_block, T_K, pred["avg_sigma_emit_eff"], self.b,
            emission_hazard_integral=cycles * pred["H_emit_per_cycle"],
        )
        self.W_emit += (pred["avg_sigma_emit_eff"] * self.b * self.f.L_pz *
                        kinetics["dN_emit"])
        self._sync_compat()

        # Same causal ordering as monotonic: state evolves, then opening samples
        # the updated field.  Re-evaluate one cycle and multiply by cycle count.
        post = self.predict_fatigue_cycle(waveform, T_K, n_phase)
        dB = cycles * post["mu_cleave_per_cycle"]
        if self.f.tau_B > 0 and dt_block > 0:
            self.B *= np.exp(-min(dt_block / self.f.tau_B, 80.0))
        self.B += dB
        self.t += dt_block
        self.K_prev = float(waveform.Kmax)
        renew = self._renew(dt_block)
        sig_peak = self.sigma_tip(float(waveform.Kmax))
        lc, lcraw, Gc = self.lambda_cleave(sig_peak, T_K)
        le, _, Ge = self.lambda_emit(sig_peak, T_K)
        mpzd = self.mpz_state.diagnostics(self.G, self.nu, self.b, self.f.r0, self.f.c_blunt)
        return {
            **renew,
            "cycles": cycles, "time_s": float(self.t),
            "Kmax_Pa_sqrt_m": float(waveform.Kmax),
            "DeltaK_Pa_sqrt_m": float(waveform.DeltaK),
            "R": float(waveform.R), "frequency_Hz": float(waveform.frequency_Hz),
            "T_K": float(T_K),
            "mu_emit": float(pred["dN_emit_per_cycle"]),
            "mu_peierls": float(kinetics.get("peierls_rate_s", 0.0) * waveform.period_s),
            "mu_taylor": float(kinetics.get("taylor_completion_rate_s", 0.0) * waveform.period_s),
            "mu_escape": float(kinetics["dN_escaped"] / max(cycles, 1.0e-300)),
            "mu_cleave_pred": float(post["mu_cleave_per_cycle"]),
            "lambda_e": float(le), "lambda_c": float(lc), "lambda_c_raw": float(lcraw),
            "G_emit_eV": float(Ge / EV_TO_J),
            "G_peierls_eV": float(kinetics.get("G_peierls_eV", 0.0)),
            "G_taylor_eV": float(kinetics.get("G_taylor_eV", 0.0)),
            "G_cleave_eff_eV": float(Gc / EV_TO_J),
            **self.cleavage_diagnostics(sig_peak, T_K),
            "sigma_tip": float(sig_peak), "sigma_back": float(self.sigma_back()),
            "r_eff": float(self.r_eff()),
            "store_per_cycle": float(kinetics["dN_trapped"] / max(cycles, 1.0e-300)),
            "mobile_per_cycle": float(kinetics["dN_emit"] / max(cycles, 1.0e-300)),
            "escape_per_cycle": float(kinetics["dN_escaped"] / max(cycles, 1.0e-300)),
            "peierls_per_cycle": float(kinetics.get("peierls_rate_s", 0.0) * waveform.period_s),
            "taylor_per_cycle": float(kinetics.get("taylor_completion_rate_s", 0.0) * waveform.period_s),
            "storage_fraction": float(kinetics["dN_trapped"] / max(kinetics["dN_emit"], 1.0e-300)),
            "dN_emit_block": float(kinetics["dN_emit"]),
            "dN_peierls_block": float(kinetics.get("peierls_rate_s", 0.0) * dt_block),
            "dN_taylor_block": float(kinetics.get("taylor_completion_rate_s", 0.0) * dt_block),
            "dN_escape_block": float(kinetics["dN_escaped"]),
            "dN_mobile_block": float(kinetics["dN_emit"]),
            "dN_store_block": float(kinetics["dN_trapped"]),
            "dN_recover_block": float(kinetics["dN_recovered"] + kinetics["dN_annihilated"]),
            "dB_block": float(dB), "B": float(self.B), "N_em": float(self.N_em),
            "r_eff_m": float(self.r_eff()), "sigma_back_Pa": float(self.sigma_back()),
            "dG_emb_eV": 0.0, "a_adv_m": float(self.a_adv), "n_adv": int(self.n_adv),
            "avg_sigma_tip_Pa": float(post["avg_sigma_tip"]),
            "max_sigma_tip_Pa": float(post["max_sigma_tip"]),
            "avg_sigma_emit_eff_Pa": float(post["avg_sigma_emit_eff"]),
            "front_state_model_code": 1.0,
            **kinetics, **mpzd,
        }

    def export_process_zone_state(self) -> dict[str, Any]:
        return self.mpz_state.state_dict()
