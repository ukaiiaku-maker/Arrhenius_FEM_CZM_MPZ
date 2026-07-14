from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    if new in text:
        return
    if old not in text:
        raise SystemExit(
            f"required patch anchor not found in {path}: {old[:80]!r}"
        )
    p.write_text(text.replace(old, new, 1))


replace_once(
    "fit_mpz_three_classes.py",
    '''        "--mpz-pair-annihilation-rate-per-count-s",\n        fs(row.mpz_pair_annihilation_rate_per_count_s),\n''',
    '''        "--mpz-pair-annihilation-rate-per-count-s",\n        fs(row.mpz_pair_annihilation_rate_per_count_s),\n        "--bulk-kinetics-model", "emission_derived_peierls_taylor_multihit",\n        "--peierls-energy-scale", fs(row.get("pt_peierls_energy_ratio", 0.005)),\n        "--peierls-entropy-scale", fs(row.get("pt_peierls_entropy_ratio", 0.005)),\n        "--taylor-energy-scale", fs(row.get("pt_taylor_energy_ratio", 0.02)),\n        "--taylor-entropy-scale", fs(row.get("pt_taylor_entropy_ratio", 0.02)),\n        "--pt-taylor-corr-rho-c", fs(row.get("pt_taylor_corr_rho_c", 1.0e14)),\n        "--pt-taylor-renewal-time-s", fs(row.get("pt_taylor_renewal_time_s", 1.0e-9)),\n        "--pt-taylor-m-exponent", fs(row.get("pt_taylor_m_exponent", 1.0)),\n        "--pt-taylor-m-scale", fs(row.get("pt_taylor_m_scale", 1.0)),\n        "--pt-taylor-m-cap", fs(row.get("pt_taylor_m_cap", float("inf"))),\n        "--pt-mobile-fraction", fs(row.get("pt_mobile_fraction", 0.01)),\n        "--pt-mobile-saturation-density-m2",\n        fs(row.get("pt_mobile_saturation_density_m2", 1.0e14)),\n''',
)

replace_once(
    "arrhenius_fracture/mpz_front_engine.py",
    '''            "mu_peierls": 0.0, "mu_taylor": 0.0,\n''',
    '''            "mu_peierls": float(kinetics.get("peierls_rate_s", 0.0) * waveform.period_s),\n            "mu_taylor": float(kinetics.get("taylor_completion_rate_s", 0.0) * waveform.period_s),\n''',
)
replace_once(
    "arrhenius_fracture/mpz_front_engine.py",
    '''            "G_peierls_eV": float(self.mpz_config.glide_barrier_eV),\n            "G_taylor_eV": float(self.mpz_config.trap_barrier_eV),\n''',
    '''            "G_peierls_eV": float(kinetics.get("G_peierls_eV", 0.0)),\n            "G_taylor_eV": float(kinetics.get("G_taylor_eV", 0.0)),\n''',
)
replace_once(
    "arrhenius_fracture/mpz_front_engine.py",
    '''            "peierls_per_cycle": 0.0, "taylor_per_cycle": 0.0,\n''',
    '''            "peierls_per_cycle": float(kinetics.get("peierls_rate_s", 0.0) * waveform.period_s),\n            "taylor_per_cycle": float(kinetics.get("taylor_completion_rate_s", 0.0) * waveform.period_s),\n''',
)
replace_once(
    "arrhenius_fracture/mpz_front_engine.py",
    '''            "dN_peierls_block": 0.0, "dN_taylor_block": 0.0,\n''',
    '''            "dN_peierls_block": float(kinetics.get("peierls_rate_s", 0.0) * dt_block),\n            "dN_taylor_block": float(kinetics.get("taylor_completion_rate_s", 0.0) * dt_block),\n''',
)
