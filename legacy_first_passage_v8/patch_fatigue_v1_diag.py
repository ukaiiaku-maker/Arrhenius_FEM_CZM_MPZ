from pathlib import Path
p=Path('/mnt/data/patch_cleave/arrhenius_fracture/fatigue_v1.py')
s=p.read_text()
# Add derivative methods in ScaledExpFloorBarrier before as_dict
marker='''    def as_dict(self) -> dict:\n'''
insert=r'''    def dG_dsigma_eV_per_GPa_numeric(self, sigma_Pa: float, T_K: float) -> float:
        """Diagnostic derivative dDeltaG/dsigma [eV/GPa]."""
        sig = abs(float(sigma_Pa))
        h = max(1.0e5, 1.0e-5 * max(sig, 1.0))
        sm = max(sig - h, 0.0)
        sp = sig + h
        Gp = float(self.deltaG_eV(sp, T_K))
        Gm = float(self.deltaG_eV(sm, T_K))
        return (Gp - Gm) / max(sp - sm, 1.0) * 1.0e9

    def vstar_b3_numeric(self, sigma_Pa: float, T_K: float, b: float = 2.74e-10) -> float:
        """Diagnostic phi v* = -dDeltaG/dsigma in b^3."""
        dG_eV_per_Pa = self.dG_dsigma_eV_per_GPa_numeric(sigma_Pa, T_K) / 1.0e9
        v_m3 = -(dG_eV_per_Pa * EV_TO_J)
        return v_m3 / max(b**3, 1.0e-300)

'''
if marker not in s: raise SystemExit('as_dict marker not found')
s=s.replace(marker,insert+marker)
# Add output diagnostics around G_emit_eV/G_cleave
old='''            "G_emit_eV": float(self.emit_barrier.deltaG_eV(pred.avg_sigma_emit_eff, T_K)),\n            "G_cleave_eff_eV": float(Gcleave_eff_diag / EV_TO_J),\n            "sigma_tip": float(pred.max_sigma_tip),\n'''
new='''            "G_emit_eV": float(self.emit_barrier.deltaG_eV(pred.avg_sigma_emit_eff, T_K)),\n            "S_emit_kB": float(self.emit_barrier.entropy_over_kB_numeric(pred.avg_sigma_emit_eff, T_K)),\n            "dGemit_dsigma_eV_per_GPa": float(self.emit_barrier.dG_dsigma_eV_per_GPa_numeric(pred.avg_sigma_emit_eff, T_K)),\n            "vstar_emit_b3": float(self.emit_barrier.vstar_b3_numeric(pred.avg_sigma_emit_eff, T_K, getattr(front, 'b', 2.74e-10))),\n            "G_peierls_eV": float(self.peierls_barrier.deltaG_eV(pred.avg_sigma_emit_eff, T_K)),\n            "S_peierls_kB": float(self.peierls_barrier.entropy_over_kB_numeric(pred.avg_sigma_emit_eff, T_K)),\n            "G_taylor_eV": float(self.taylor_barrier.deltaG_eV(pred.avg_sigma_emit_eff, T_K)),\n            "S_taylor_kB": float(self.taylor_barrier.entropy_over_kB_numeric(pred.avg_sigma_emit_eff, T_K)),\n            "G_cleave_eff_eV": float(Gcleave_eff_diag / EV_TO_J),\n            **front.cleavage_diagnostics(pred.max_sigma_tip, T_K),\n            "sigma_tip": float(pred.max_sigma_tip),\n'''
if old not in s: raise SystemExit('out diag marker not found')
s=s.replace(old,new)
p.write_text(s)
