from pathlib import Path
p=Path('/mnt/data/patch_cleave/arrhenius_fracture/sharp_front.py')
s=p.read_text()
# Insert helper to apply cleavage barrier args before default_cleavage_barrier
marker='''def default_cleavage_barrier() -> FractureBarrier:\n'''
helper=r'''
def apply_cleavage_barrier_args(cb: FractureBarrier, args) -> FractureBarrier:
    """Apply CLI/namespace cleavage-barrier overrides.

    The code now supports both the legacy H-TS-sigma*v cleavage barrier
    (``classic``) and a direct EXP-floor free-energy surface
    (``exp_floor``).  In exp_floor mode entropy and activation volume are
    diagnostics from derivatives of DeltaG*(sigma,T), not independent terms.
    """
    # Backward-compatible classic overrides.
    if getattr(args, 'cleave_H0_eV', None) is not None:
        cb.H0_eV = float(args.cleave_H0_eV)
    if getattr(args, 'cleave_S0_kB', None) is not None:
        cb.S0_neg_kB = float(args.cleave_S0_kB)
    if getattr(args, 'cleave_sigma0_S', None) is not None:
        cb.sigma0_S_GPa = float(args.cleave_sigma0_S)
    if getattr(args, 'cleave_S_sigma_max_kB', None) is not None:
        cb.S_sigma_max_kB = float(args.cleave_S_sigma_max_kB)
    if getattr(args, 'cleave_entropy_form', None) is not None:
        cb.entropy_stress_form = str(args.cleave_entropy_form)
    # Shared entropy-form override used by the old sweep code.
    if getattr(args, 'entropy_form', None) is not None:
        cb.entropy_stress_form = str(args.entropy_form)
        cb.use_negative_entropy = True
    if getattr(args, 'entropy_gate_power', None) is not None:
        cb.entropy_gate_power = float(args.entropy_gate_power)

    kind = getattr(args, 'cleave_barrier_kind', None)
    if kind is not None:
        cb.barrier_kind = str(kind)
    # EXP-floor controls.  These are no-ops unless barrier_kind == exp_floor,
    # but we always populate them so run_args.json fully records the surface.
    mapping = {
        'cleave_G00_eV': 'ef_G00_eV',
        'cleave_gT_eV_per_K': 'ef_gT_eV_per_K',
        'cleave_sigc0_GPa': 'ef_sigc0_Pa',
        'cleave_sT_GPa_per_K': 'ef_sT_Pa_per_K',
        'cleave_exp_a': 'ef_a',
        'cleave_exp_n': 'ef_n',
        'cleave_floor_frac': 'ef_floor_frac',
        'cleave_floor_min_eV': 'ef_floor_min_eV',
        'cleave_floor_max_frac': 'ef_floor_max_frac',
        'cleave_Tref_K': 'ef_Tref_K',
        'cleave_exp_T_mode': 'ef_T_mode',
        'cleave_mu_dlnmu_dT_per_K': 'ef_mu_dlnmu_dT_per_K',
        'cleave_G0_mu_power': 'ef_G0_mu_power',
        'cleave_sigc_mu_power': 'ef_sigc_mu_power',
        'cleave_S_hs_kB': 'ef_S_hs_kB',
        'cleave_sigma_S_GPa': 'ef_sigma_S_GPa',
        'cleave_S_hs_power': 'ef_S_hs_power',
        'cleave_S_hs_dT_per_K': 'ef_S_hs_dT_per_K',
        'cleave_S_hs_Tref_K': 'ef_S_hs_Tref_K',
    }
    for src, dst in mapping.items():
        if getattr(args, src, None) is None:
            continue
        val = getattr(args, src)
        if src in ('cleave_sigc0_GPa', 'cleave_sT_GPa_per_K', 'cleave_sigma_S_GPa'):
            val = float(val) * 1.0e9
        setattr(cb, dst, val)
    if getattr(args, 'cleave_monotone_stress', None) is not None:
        cb.monotone_stress = bool(args.cleave_monotone_stress)
    return cb

'''
if marker not in s: raise SystemExit('marker default not found')
s=s.replace(marker, helper+marker)
# Replace manual cb overrides in build_engine with apply helper. Find block.
old='''    cb = default_cleavage_barrier()\n    if args.cleave_H0_eV is not None: cb.H0_eV = args.cleave_H0_eV\n    if args.cleave_S0_kB is not None: cb.S0_neg_kB = args.cleave_S0_kB\n    if getattr(args, 'cleave_sigma0_S', None) is not None:\n        cb.sigma0_S_GPa = args.cleave_sigma0_S\n    eb = default_emission_barrier(mat.b)\n'''
new='''    cb = apply_cleavage_barrier_args(default_cleavage_barrier(), args)\n    eb = default_emission_barrier(mat.b)\n'''
if old not in s: raise SystemExit('build_engine cb block not found')
s=s.replace(old,new)
# Remove later duplicate cleave_S_sigma override? It now exists; but okay setting twice? The later still exists. Remove to avoid inconsistency? It comes after entropy_form with cleave_S_sigma. We'll leave harmless.
# Add diagnostics in lambda_cleave? Add new method after lambda_cleave.
marker='''    def predict_clock_increment(self, K, T, dt):\n'''
insert=r'''    def cleavage_diagnostics(self, sig_tip, T):
        """Return raw/effective cleavage free-energy diagnostics at sig_tip."""
        sig_eff = max(float(sig_tip) - self.f.chi_shield * self.sigma_back(), 0.0)
        d = self.cb.diagnostics(np.array([sig_eff]), T, self.b)
        Gstar = float(d['G_eV'][0] * EV_TO_J)
        dGe = min(self.dG_emb(), self.f.emb_sat_frac * Gstar)
        Geff = max(Gstar - dGe, 0.0)
        return {
            'sigma_cleave_eff_Pa': float(sig_eff),
            'G_cleave_raw_eV': float(Gstar / EV_TO_J),
            'G_cleave_eff_eV': float(Geff / EV_TO_J),
            'S_cleave_kB': float(d['S_kB'][0]),
            'dGcleave_dsigma_eV_per_GPa': float(d['dG_dsigma_eV_per_GPa'][0]),
            'vstar_cleave_b3': float(d['vstar_b3'][0]),
            'cleave_barrier_kind_code': 1.0 if str(getattr(self.cb, 'barrier_kind', 'classic')) == 'exp_floor' else 0.0,
        }

'''
if marker not in s: raise SystemExit('predict marker not found')
s=s.replace(marker,insert+marker)
# In fatigue branch after Gcleave_eff diag in fatigue_v1 handles. For sharp FrontEngine.step add diagnostics to info near G_cleave_eff_eV.
old="""            'dG_emb_eV': self.dG_emb() / EV_TO_J, 'G_cleave_eff_eV': Gc_eff / EV_TO_J,\n"""
new="""            'dG_emb_eV': self.dG_emb() / EV_TO_J, 'G_cleave_eff_eV': Gc_eff / EV_TO_J,\n            **self.cleavage_diagnostics(sig2, T),\n"""
if old not in s: raise SystemExit('step info G_cleave line not found')
s=s.replace(old,new)
# In 2D rows append add diagnostics after cyclic_plastic_work_acc
old="""                         pz_store_total, pz_mobile_total, pz_escape_total, pz_emit_total,\n                         int(cyclic_mechanics_updates), cyclic_plastic_work_acc))\n"""
new="""                         pz_store_total, pz_mobile_total, pz_escape_total, pz_emit_total,\n                         int(cyclic_mechanics_updates), cyclic_plastic_work_acc,\n                         float(info.get('G_cleave_raw_eV', info.get('G_cleave_eff_eV', 0.0))),\n                         float(info.get('G_cleave_eff_eV', 0.0)),\n                         float(info.get('S_cleave_kB', 0.0)),\n                         float(info.get('dGcleave_dsigma_eV_per_GPa', 0.0)),\n                         float(info.get('vstar_cleave_b3', 0.0)),\n                         float(info.get('sigma_cleave_eff_Pa', info.get('sigma_tip', 0.0))),\n                         float(info.get('cleave_barrier_kind_code', 0.0))))\n"""
if old not in s: raise SystemExit('rows append marker not found')
s=s.replace(old,new)
# Header add columns at end
old="""                          'cyclic_mechanics_updates,cyclic_plastic_work_J',\n"""
new="""                          'cyclic_mechanics_updates,cyclic_plastic_work_J,'\n                          'G_cleave_raw_eV,G_cleave_eff_eV,S_cleave_kB,'\n                          'dGcleave_dsigma_eV_per_GPa,vstar_cleave_b3,'\n                          'sigma_cleave_eff_Pa,cleave_barrier_kind_code',\n"""
if old not in s: raise SystemExit('header marker not found')
s=s.replace(old,new)
# Add parser args after cleave-S-sigma
old="""    p.add_argument('--cleave-S-sigma-max-kB', type=float, default=None,\n                   dest='cleave_S_sigma_max_kB',\n                   help='Cleavage Schoeck stress-term magnitude [kB] (physical form)')\n\n    # sweep mode: each axis is a list; the Cartesian product is swept, and each\n"""
new="""    p.add_argument('--cleave-S-sigma-max-kB', type=float, default=None,\n                   dest='cleave_S_sigma_max_kB',\n                   help='Cleavage Schoeck stress-term magnitude [kB] (physical form)')\n    p.add_argument('--cleave-barrier-kind', choices=['classic', 'exp_floor'], default=None,\n                   dest='cleave_barrier_kind',\n                   help='Cleavage free-energy surface. classic uses H-TS-sigma*v; exp_floor uses DeltaG*(sigma,T) directly.')\n    p.add_argument('--cleave-G00-eV', type=float, default=None, dest='cleave_G00_eV')\n    p.add_argument('--cleave-gT-eV-per-K', type=float, default=None, dest='cleave_gT_eV_per_K')\n    p.add_argument('--cleave-sigc0-GPa', type=float, default=None, dest='cleave_sigc0_GPa')\n    p.add_argument('--cleave-sT-GPa-per-K', type=float, default=None, dest='cleave_sT_GPa_per_K')\n    p.add_argument('--cleave-exp-a', type=float, default=None, dest='cleave_exp_a')\n    p.add_argument('--cleave-exp-n', type=float, default=None, dest='cleave_exp_n')\n    p.add_argument('--cleave-floor-frac', type=float, default=None, dest='cleave_floor_frac')\n    p.add_argument('--cleave-floor-min-eV', type=float, default=None, dest='cleave_floor_min_eV')\n    p.add_argument('--cleave-floor-max-frac', type=float, default=None, dest='cleave_floor_max_frac')\n    p.add_argument('--cleave-Tref-K', type=float, default=None, dest='cleave_Tref_K')\n    p.add_argument('--cleave-exp-T-mode', choices=['linear', 'mu_scale'], default=None, dest='cleave_exp_T_mode',\n                   help='Temperature law for exp_floor cleavage: linear or shear-modulus-like mu_scale.')\n    p.add_argument('--cleave-mu-dlnmu-dT-per-K', type=float, default=None, dest='cleave_mu_dlnmu_dT_per_K')\n    p.add_argument('--cleave-G0-mu-power', type=float, default=None, dest='cleave_G0_mu_power')\n    p.add_argument('--cleave-sigc-mu-power', type=float, default=None, dest='cleave_sigc_mu_power')\n    p.add_argument('--cleave-S-hs-kB', type=float, default=None, dest='cleave_S_hs_kB',\n                   help='Optional high-stress entropy shift for exp_floor [kB].')\n    p.add_argument('--cleave-sigma-S-GPa', type=float, default=None, dest='cleave_sigma_S_GPa')\n    p.add_argument('--cleave-S-hs-power', type=float, default=None, dest='cleave_S_hs_power')\n    p.add_argument('--cleave-S-hs-dT-per-K', type=float, default=None, dest='cleave_S_hs_dT_per_K')\n    p.add_argument('--cleave-S-hs-Tref-K', type=float, default=None, dest='cleave_S_hs_Tref_K')\n    p.add_argument('--cleave-monotone-stress', action=argparse.BooleanOptionalAction, default=None, dest='cleave_monotone_stress',\n                   help='Apply monotone-stress envelope to classic cleavage barrier diagnostics/rates. Default uses barrier preset.')\n\n    # sweep mode: each axis is a list; the Cartesian product is swept, and each\n"""
if old not in s: raise SystemExit('parser cleave insertion not found')
s=s.replace(old,new)
p.write_text(s)
