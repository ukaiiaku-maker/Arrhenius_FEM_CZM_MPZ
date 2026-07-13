"""Shared command-line serialization for v9 moving-PZ parameter tables."""
from __future__ import annotations

import math
import sys


def fs(x):
    x=float(x)
    return 'inf' if math.isinf(x) else f'{x:.16g}'


def moving_pz_cli(row):
    return [
        '--front-state-model','moving_pz','--sigma-cap-GPa','0',
        '--r-pz',fs(row.r_pz_m),'--c-blunt',fs(row.c_blunt),'--L-pz',fs(row.mpz_length_m),
        '--mpz-length-m',fs(row.mpz_length_m),'--mpz-n-bins',str(int(row.mpz_n_bins)),
        '--mpz-n-systems',str(int(row.mpz_n_systems)),
        '--mpz-source-sites-per-system',fs(row.mpz_source_sites_per_system),
        '--mpz-source-recovery-rate-s',fs(row.mpz_source_recovery_rate_s),
        '--mpz-source-refresh-length-m',fs(row.mpz_source_refresh_length_m),
        '--mpz-shielding-factors',str(row.mpz_shielding_factors),
        '--mpz-glide-barrier-eV',fs(row.mpz_glide_barrier_eV),
        '--mpz-glide-activation-volume-b3',fs(row.mpz_glide_activation_volume_b3),
        '--mpz-trap-barrier-eV',fs(row.mpz_trap_barrier_eV),
        '--mpz-detrap-barrier-eV',fs(row.mpz_detrap_barrier_eV),
        '--mpz-retained-recovery-barrier-eV',fs(row.mpz_retained_recovery_barrier_eV),
        '--mpz-pair-annihilation-rate-per-count-s',fs(row.mpz_pair_annihilation_rate_per_count_s),
        '--cleave-barrier-kind','exp_floor','--cleave-exp-T-mode','linear',
        '--cleave-G00-eV',fs(row.cleave_G00_eV),'--cleave-gT-eV-per-K',fs(row.cleave_gT_eV_per_K),
        '--cleave-sigc0-GPa',fs(row.cleave_sigc0_GPa),'--cleave-sT-GPa-per-K',fs(row.cleave_sT_GPa_per_K),
        '--cleave-exp-a',fs(row.cleave_exp_a),'--cleave-exp-n',fs(row.cleave_exp_n),
        '--cleave-floor-frac',fs(row.cleave_floor_frac),'--cleave-S-hs-kB',fs(row.cleave_S_hs_kB),
        '--cleave-sigma-S-GPa','6','--cleave-S-hs-power','2',
        '--emit-barrier-kind','exp_floor','--emit-G00-eV',fs(row.emit_G00_eV),
        '--emit-gT-eV-per-K',fs(row.emit_gT_eV_per_K),'--emit-sigc0-GPa',fs(row.emit_sigc0_GPa),
        '--emit-sT-GPa-per-K',fs(row.emit_sT_GPa_per_K),'--emit-exp-a',fs(row.emit_exp_a),
        '--emit-exp-n',fs(row.emit_exp_n),'--emit-floor-frac',fs(row.emit_floor_frac),
    ]


def check_parameter_status(table, source, require_fitted=False):
    """Warn or fail when a runner is given uncalibrated MPZ rows."""
    if 'status' not in table.columns:
        msg=f"parameter table {source} has no status column"
        if require_fitted:
            raise SystemExit(msg)
        print(f"WARNING: {msg}", file=sys.stderr)
        return
    statuses=sorted(set(str(x) for x in table['status'].fillna('UNKNOWN')))
    unfitted=[x for x in statuses if not x.startswith('FITTED_MPZ_V9')]
    if unfitted:
        msg=(f"parameter table {source} contains non-production rows: {unfitted}. "
             "These are suitable for smoke tests only.")
        if require_fitted:
            raise SystemExit(msg)
        print(f"WARNING: {msg}", file=sys.stderr)
