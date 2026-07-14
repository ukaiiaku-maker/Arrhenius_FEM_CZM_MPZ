from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    if new in text:
        return
    if old not in text:
        raise SystemExit(f"patch anchor not found in {path}: {old[:120]!r}")
    p.write_text(text.replace(old, new, 1))


path = "arrhenius_fracture/emission_derived_plasticity.py"

replace_once(
    path,
    '''    def peierls_rate(
        self, sigma_eq_Pa: np.ndarray | float, T_K: float
    ) -> np.ndarray:
        tau = max(
            float(self.cfg.peierls_stress_fraction), 0.0
        ) * np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
        G = self.barrier_eV("peierls", tau, T_K)
        return self._arrhenius_rate(
            G, T_K, self.cfg.peierls.rate_prefactor_s
        )
''',
    '''    def raw_zero_stress_barrier_eV(
        self, mechanism: str, T_K: float
    ) -> float:
        """Unclamped zero-stress free energy used for parameter validation."""
        key = str(mechanism).lower()
        if key in {"peierls", "p"}:
            scale = self.cfg.peierls
        elif key in {"taylor", "t"}:
            scale = self.cfg.taylor
        elif key in {"emission", "emit", "e"}:
            scale = MechanismScale(1.0, 1.0, 1.0, 1.0)
        else:
            raise ValueError(f"unknown mechanism: {mechanism}")
        dT = float(T_K) - float(self.cfg.parent.Tref_K)
        return (
            float(scale.energy_ratio) * float(self.cfg.parent.G00_eV)
            + float(scale.entropy_ratio)
            * float(self.cfg.parent.gT_eV_per_K) * dT
        )

    def peierls_rates(
        self, sigma_eq_Pa: np.ndarray | float, T_K: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Net, forward, and reverse Peierls rates.

        The reverse reference uses the zero-stress EXP-floor barrier, matching
        the signed detailed-balance Peierls law used in the prior DDD model.
        Consequently forward and reverse rates are identical at zero stress and
        the net plastic drift is exactly zero.
        """
        tau = max(
            float(self.cfg.peierls_stress_fraction), 0.0
        ) * np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
        G_forward = self.barrier_eV("peierls", tau, T_K)
        G_reverse = self.barrier_eV("peierls", 0.0, T_K)
        forward = self._arrhenius_rate(
            G_forward, T_K, self.cfg.peierls.rate_prefactor_s
        )
        reverse = self._arrhenius_rate(
            G_reverse, T_K, self.cfg.peierls.rate_prefactor_s
        )
        net = np.maximum(forward - reverse, 0.0)
        return net, forward, np.broadcast_to(reverse, np.shape(forward))

    def peierls_rate(
        self, sigma_eq_Pa: np.ndarray | float, T_K: float
    ) -> np.ndarray:
        return self.peierls_rates(sigma_eq_Pa, T_K)[0]
''',
)

replace_once(
    path,
    '''    def taylor_single_hit_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> np.ndarray:
        sigma_local = self.taylor_local_stress(
            sigma_eq_Pa, rho_forest_m2, b_m
        )
        G = self.barrier_eV("taylor", sigma_local, T_K)
        return self._arrhenius_rate(
            G, T_K, self.cfg.taylor.rate_prefactor_s
        )

    def taylor_completion_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        h1 = self.taylor_single_hit_rate(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        m = self.cfg.correlated_taylor.hit_order(rho_forest_m2)
        tc = max(
            float(self.cfg.correlated_taylor.renewal_time_s), 1.0e-30
        )
        exposure = np.minimum(np.maximum(h1 * tc, 0.0), 1.0e12)
        rate = gammainc(m, exposure) / tc
        return np.maximum(rate, 0.0), m, h1
''',
    '''    def taylor_single_hit_rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Net, forward, and zero-stress reverse Taylor hit rates."""
        sigma_local = self.taylor_local_stress(
            sigma_eq_Pa, rho_forest_m2, b_m
        )
        G_forward = self.barrier_eV("taylor", sigma_local, T_K)
        G_reverse = self.barrier_eV("taylor", 0.0, T_K)
        forward = self._arrhenius_rate(
            G_forward, T_K, self.cfg.taylor.rate_prefactor_s
        )
        reverse = self._arrhenius_rate(
            G_reverse, T_K, self.cfg.taylor.rate_prefactor_s
        )
        reverse = np.broadcast_to(reverse, np.shape(forward))
        return np.maximum(forward - reverse, 0.0), forward, reverse

    def taylor_single_hit_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> np.ndarray:
        return self.taylor_single_hit_rates(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )[0]

    def taylor_completion_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        single_net, single_forward, single_reverse = (
            self.taylor_single_hit_rates(
                sigma_eq_Pa, rho_forest_m2, T_K, b_m
            )
        )
        m = self.cfg.correlated_taylor.hit_order(rho_forest_m2)
        tc = max(
            float(self.cfg.correlated_taylor.renewal_time_s), 1.0e-30
        )
        exposure_forward = np.minimum(
            np.maximum(single_forward * tc, 0.0), 1.0e12
        )
        exposure_reverse = np.minimum(
            np.maximum(single_reverse * tc, 0.0), 1.0e12
        )
        completion_forward = gammainc(m, exposure_forward) / tc
        completion_reverse = gammainc(m, exposure_reverse) / tc
        completion_net = np.maximum(
            completion_forward - completion_reverse, 0.0
        )
        return (
            completion_net, m, single_net, single_forward, single_reverse,
            completion_forward, completion_reverse,
        )
''',
)

replace_once(
    path,
    '''        p = self.peierls_rate(sigma_eq_Pa, T_K)
        t, m, t1 = self.taylor_completion_rate(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        series = self.series_rate(p, t)
''',
    '''        p, p_forward, p_reverse = self.peierls_rates(
            sigma_eq_Pa, T_K
        )
        (
            t, m, t1, t1_forward, t1_reverse,
            t_forward, t_reverse,
        ) = self.taylor_completion_rate(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        series = self.series_rate(p, t)
''',
)

replace_once(
    path,
    '''            "peierls_rate_s": np.asarray(p, dtype=float),
            "taylor_single_hit_rate_s": np.asarray(t1, dtype=float),
            "taylor_completion_rate_s": np.asarray(t, dtype=float),
''',
    '''            "peierls_rate_s": np.asarray(p, dtype=float),
            "peierls_forward_rate_s": np.asarray(p_forward, dtype=float),
            "peierls_reverse_rate_s": np.asarray(p_reverse, dtype=float),
            "taylor_single_hit_rate_s": np.asarray(t1, dtype=float),
            "taylor_single_hit_forward_rate_s": np.asarray(t1_forward, dtype=float),
            "taylor_single_hit_reverse_rate_s": np.asarray(t1_reverse, dtype=float),
            "taylor_completion_rate_s": np.asarray(t, dtype=float),
            "taylor_completion_forward_rate_s": np.asarray(t_forward, dtype=float),
            "taylor_completion_reverse_rate_s": np.asarray(t_reverse, dtype=float),
''',
)

sp = "search_mpz_peierls_taylor_parameters.py"
replace_once(
    sp,
    '"pt_entropy_multiplier": log_scale(u[:, 2], 0.25, 128.0),',
    '"pt_entropy_multiplier": log_scale(u[:, 2], 0.25, 8.0),',
)
replace_once(
    sp,
    '''    curves = []
    slopes = []
    drops = []
    resolved = True
''',
    '''    curves = []
    slopes = []
    drops = []
    resolved = True
    raw_peierls = np.asarray([
        model.raw_zero_stress_barrier_eV("peierls", T)
        for T in temperatures
    ], dtype=float)
    raw_taylor = np.asarray([
        model.raw_zero_stress_barrier_eV("taylor", T)
        for T in temperatures
    ], dtype=float)
    min_raw_barrier = float(min(raw_peierls.min(), raw_taylor.min()))
    barrier_valid = bool(min_raw_barrier > 1.0e-8)
    zero_rate_max = 0.0
    for T in temperatures:
        zr = model.rates(0.0, rho, T, b)["equivalent_plastic_rate_s"]
        zero_rate_max = max(zero_rate_max, float(np.max(np.abs(zr))))
    detailed_balance_valid = bool(zero_rate_max <= 1.0e-20)
''',
)
replace_once(
    sp,
    '''    admissible = (
        resolved
        and monotonic
        and np.isfinite(smax)
        and smax <= max_stress_GPa
    )
''',
    '''    admissible = (
        resolved
        and monotonic
        and barrier_valid
        and detailed_balance_valid
        and np.isfinite(smax)
        and smax <= max_stress_GPa
    )
''',
)
replace_once(
    sp,
    '''    if strength_window:
        status = "strict_strength_window"
    elif admissible:
        status = "monotonic_topology_only"
    elif not resolved:
        status = "unresolved_at_sigma_limit"
    elif not monotonic:
        status = "high_density_downturn"
    else:
        status = "stress_above_limit"
''',
    '''    if strength_window:
        status = "strict_strength_window"
    elif admissible:
        status = "monotonic_topology_only"
    elif not barrier_valid:
        status = "invalid_scaled_zero_stress_barrier"
    elif not detailed_balance_valid:
        status = "zero_stress_ratchet"
    elif not resolved:
        status = "unresolved_at_sigma_limit"
    elif not monotonic:
        status = "high_density_downturn"
    else:
        status = "stress_above_limit"
''',
)
replace_once(
    sp,
    '''        "resolved": bool(resolved),
        "monotonic": bool(monotonic),
        "accepted": bool(admissible),
''',
    '''        "resolved": bool(resolved),
        "monotonic": bool(monotonic),
        "barrier_valid": bool(barrier_valid),
        "detailed_balance_valid": bool(detailed_balance_valid),
        "min_raw_scaled_G0_eV": min_raw_barrier,
        "zero_stress_rate_max_s": zero_rate_max,
        "accepted": bool(admissible),
''',
)
replace_once(
    sp,
    'ap.add_argument("--max-stress-GPa", type=float, default=80.0)',
    'ap.add_argument("--max-stress-GPa", type=float, default=40.0)',
)
replace_once(
    sp,
    '''            stress_penalty = abs(math.log(sigma_ref / 2.0))
''',
    '''            stress_penalty = (
                abs(math.log(sigma_ref / 2.0))
                + max(float(metrics["sigma_max_GPa"]) - 20.0, 0.0) / 10.0
            )
''',
)

jp = "prepare_mpz_v9_3_pt_input.py"
replace_once(
    jp,
    '''import argparse
from pathlib import Path

import pandas as pd
''',
    '''import argparse
from pathlib import Path
import warnings

import pandas as pd
''',
)
replace_once(
    jp,
    '''    joined = atlas.merge(
        materials,
        on="candidate_id",
        how="left",
        suffixes=("", "__material"),
        validate="many_to_one",
        indicator=True,
    )
''',
    '''    with warnings.catch_warnings():
        warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
        joined = atlas.merge(
            materials,
            on="candidate_id",
            how="left",
            suffixes=("", "__material"),
            validate="many_to_one",
            indicator=True,
        )
''',
)

rp = "run_mpz_v9_3_peierls_taylor_search.sh"
replace_once(
    rp,
    '''RHO_POINTS="${RHO_POINTS:-65}"
SEED="${SEED:-93017}"
''',
    '''RHO_POINTS="${RHO_POINTS:-65}"
MAX_STRESS_GPA="${MAX_STRESS_GPA:-40}"
SEED="${SEED:-93017}"
''',
)
replace_once(
    rp,
    '''echo "[$(stamp)] rho=$RHO_MIN..$RHO_MAX points=$RHO_POINTS"
''',
    '''echo "[$(stamp)] rho=$RHO_MIN..$RHO_MAX points=$RHO_POINTS max_stress_GPa=$MAX_STRESS_GPA"
''',
)
replace_once(
    rp,
    '''  --rho-points "$RHO_POINTS" \\
  --seed "$SEED" \\
''',
    '''  --rho-points "$RHO_POINTS" \\
  --max-stress-GPa "$MAX_STRESS_GPA" \\
  --seed "$SEED" \\
''',
)

tp = Path("tests/test_emission_derived_peierls_taylor.py")
text = tp.read_text()
extra = '''


def test_signed_detailed_balance_gives_zero_net_rate_at_zero_stress():
    model = EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=1.0e11, renewal_time_s=1.0e-10, m_cap=22.0
            )
        )
    )
    rho = np.logspace(12, 17, 20)
    zero = model.rates(0.0, rho, 700.0, 2.74e-10)
    assert np.all(zero["peierls_rate_s"] == 0.0)
    assert np.all(zero["taylor_single_hit_rate_s"] == 0.0)
    assert np.all(zero["taylor_completion_rate_s"] == 0.0)
    assert np.all(zero["series_rate_s"] == 0.0)
    assert np.all(zero["equivalent_plastic_rate_s"] == 0.0)

    driven = model.rates(2.0e9, rho, 700.0, 2.74e-10)
    assert np.all(
        driven["peierls_forward_rate_s"]
        >= driven["peierls_reverse_rate_s"]
    )
    assert np.all(
        driven["taylor_completion_forward_rate_s"]
        >= driven["taylor_completion_reverse_rate_s"]
    )
    assert np.any(driven["series_rate_s"] > 0.0)
'''
if "test_signed_detailed_balance_gives_zero_net_rate_at_zero_stress" not in text:
    tp.write_text(text.rstrip() + extra + "\n")

replace_once("pyproject.toml", 'version = "0.9.3"', 'version = "0.9.4"')
replace_once(
    "arrhenius_fracture/__init__.py",
    "__version__ = '0.9.3'",
    "__version__ = '0.9.4'",
)
