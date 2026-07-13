from pathlib import Path
p=Path('/mnt/data/patch_cleave/arrhenius_fracture/config.py')
s=p.read_text()
# Add dataclass fields after ef_floor_max_frac
old="""    ef_floor_max_frac: float = 0.95\n    # --- high-stress entropy crossover (fatigue-paper hypothesis) ----------"""
new="""    ef_floor_max_frac: float = 0.95\n    # Temperature mode for the EXP-floor free-energy surface.\n    #   linear:   G0=G00+gT*(T-Tref), sigc=sigc0+sT*(T-Tref).\n    #   mu_scale: G0=G00*[mu(T)/mu(Tref)]^pG and\n    #             sigc=sigc0*[mu(T)/mu(Tref)]^psig, with a simple\n    #             local shear-modulus proxy mu/mu_ref = 1+dlnmu_dT*(T-Tref).\n    # The mu_scale option is intended for cleavage sweeps where we want the\n    # temperature dependence to resemble elastic modulus softening instead of\n    # importing a nanopillar nucleation entropy slope directly.\n    ef_T_mode: str = 'linear'            # 'linear' | 'mu_scale'\n    ef_mu_dlnmu_dT_per_K: float = -1.5e-4\n    ef_G0_mu_power: float = 1.0\n    ef_sigc_mu_power: float = 1.0\n    # --- high-stress entropy crossover (fatigue-paper hypothesis) ----------"""
if old not in s:
    raise SystemExit('field insertion marker not found')
s=s.replace(old,new)
# Replace G0/sigc section in _exp_floor
old="""        G0 = self.ef_G00_eV + self.ef_gT_eV_per_K * (T - self.ef_Tref_K)\n        sigc = self.ef_sigc0_Pa + self.ef_sT_Pa_per_K * (T - self.ef_Tref_K)\n        G0 = max(G0, 1e-9); sigc = max(sigc, 1.0)\n        a, n = float(self.ef_a), float(self.ef_n)\n"""
new="""        Tmode = str(getattr(self, 'ef_T_mode', 'linear')).lower()\n        if Tmode in ('mu', 'mu_scale', 'modulus', 'shear_modulus'):\n            # Minimal shear-modulus-like temperature law for cleavage sweeps.\n            # dlnmu_dT < 0 gives softening with T.  Clamp away from zero so\n            # exploratory high-T sweeps remain numerically well posed.\n            dln = float(getattr(self, 'ef_mu_dlnmu_dT_per_K', -1.5e-4))\n            mu_ratio = max(0.05, 1.0 + dln * (T - self.ef_Tref_K))\n            pG = float(getattr(self, 'ef_G0_mu_power', 1.0))\n            ps = float(getattr(self, 'ef_sigc_mu_power', 1.0))\n            G0 = self.ef_G00_eV * (mu_ratio ** pG)\n            sigc = self.ef_sigc0_Pa * (mu_ratio ** ps)\n            dmu_ratio_dT = dln\n            dG0_dT = self.ef_G00_eV * pG * (mu_ratio ** (pG - 1.0)) * dmu_ratio_dT\n            dsigc_dT = self.ef_sigc0_Pa * ps * (mu_ratio ** (ps - 1.0)) * dmu_ratio_dT\n        else:\n            G0 = self.ef_G00_eV + self.ef_gT_eV_per_K * (T - self.ef_Tref_K)\n            sigc = self.ef_sigc0_Pa + self.ef_sT_Pa_per_K * (T - self.ef_Tref_K)\n            dG0_dT = self.ef_gT_eV_per_K\n            dsigc_dT = self.ef_sT_Pa_per_K\n        G0 = max(G0, 1e-9); sigc = max(sigc, 1.0)\n        a, n = float(self.ef_a), float(self.ef_n)\n"""
if old not in s:
    raise SystemExit('G0 block marker not found')
s=s.replace(old,new)
# Replace dExp and derivatives uses
old="""        # d(expTerm)/dT through sigc(T):  = expT*(a*n)*x^n*sigc/sigc^2 * sT\n        dExp_dT = expT * (a * n) * xn * (sigc / sigc**2) * self.ef_sT_Pa_per_K\n"""
new="""        # d(expTerm)/dT through sigc(T): exp[-a(s/sigc)^n]\n        # derivative = expT*(a*n)*x^n*(1/sigc)*dsigc/dT.\n        dExp_dT = expT * (a * n) * xn * (1.0 / sigc) * dsigc_dT\n"""
if old not in s:
    raise SystemExit('dExp marker not found')
s=s.replace(old,new)
old="""            dFloor = self.ef_floor_frac * self.ef_gT_eV_per_K\n        elif raw_floor > self.ef_floor_max_frac * G0:\n            dFloor = self.ef_floor_max_frac * self.ef_gT_eV_per_K\n"""
new="""            dFloor = self.ef_floor_frac * dG0_dT\n        elif raw_floor > self.ef_floor_max_frac * G0:\n            dFloor = self.ef_floor_max_frac * dG0_dT\n"""
if old not in s:
    raise SystemExit('dFloor marker not found')
s=s.replace(old,new)
old="""        dAmp = self.ef_gT_eV_per_K - dFloor\n"""
new="""        dAmp = dG0_dT - dFloor\n"""
if old not in s:
    raise SystemExit('dAmp marker not found')
s=s.replace(old,new)
# Add diagnostic methods before G_barrier
marker="""    def G_barrier(self, sigma: np.ndarray, T: float = 0.0, b: float = 2.74e-10) -> np.ndarray:\n"""
insert="""    def dG_dsigma_numeric(self, sigma: np.ndarray | float, T: float = 0.0,\n                          b: float = 2.74e-10) -> np.ndarray:\n        \"\"\"Numerical derivative dG*/dsigma [J/Pa] at fixed T.\n\n        This is used only for diagnostics: v* = -dG*/dsigma and the\n        fatigue-paper stationarity audit dG*/dsigma ~ 0.  It intentionally\n        differentiates the same free-energy surface used by the rate law,\n        including EXP-floor and any monotone envelope.\n        \"\"\"\n        sig = np.abs(np.asarray(sigma, dtype=float))\n        out = np.empty_like(sig, dtype=float)\n        for idx, val in np.ndenumerate(sig):\n            h = max(1.0e5, 1.0e-5 * max(float(val), 1.0))\n            sm = max(float(val) - h, 0.0)\n            sp = float(val) + h\n            Gp = float(self.G_barrier(np.array([sp]), T, b)[0])\n            Gm = float(self.G_barrier(np.array([sm]), T, b)[0])\n            denom = max(sp - sm, 1.0)\n            out[idx] = (Gp - Gm) / denom\n        return out\n\n    def diagnostics(self, sigma: np.ndarray | float, T: float = 0.0,\n                    b: float = 2.74e-10) -> dict:\n        \"\"\"Free-energy barrier diagnostics at fixed stress and temperature.\"\"\"\n        sig = np.asarray(sigma, dtype=float)\n        G = self.G_barrier(sig, T, b)\n        S = self.S(sig, T)\n        dGds = self.dG_dsigma_numeric(sig, T, b)\n        vstar = -dGds\n        return {\n            'G_eV': G / EV_TO_J,\n            'S_kB': S / KB,\n            'dG_dsigma_eV_per_GPa': dGds / EV_TO_J * 1.0e9,\n            'vstar_b3': vstar / max(b**3, 1.0e-300),\n        }\n\n"""
if marker not in s:
    raise SystemExit('G_barrier marker not found')
s=s.replace(marker, insert+marker)
p.write_text(s)
