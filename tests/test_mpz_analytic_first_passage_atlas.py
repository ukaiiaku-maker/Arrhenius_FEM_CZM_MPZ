from __future__ import annotations

import numpy as np
import pandas as pd

import build_mpz_analytic_first_passage_atlas as atlas


def shape_row() -> pd.Series:
    return pd.Series({
        "cleave_sT_GPa_per_K": 0.0,
        "cleave_exp_a": 1.2,
        "cleave_exp_n": 1.5,
        "cleave_floor_frac": 0.01,
        "emit_sT_GPa_per_K": 0.0,
        "emit_exp_a": 0.35,
        "emit_exp_n": 1.1,
        "emit_floor_frac": 0.04,
        "r_pz_m": 1.0e-6,
        "mpz_n_systems": 2,
        "mpz_source_sites_per_system": 5.0,
    })


def candidate() -> pd.DataFrame:
    return pd.DataFrame([{
        "cleave_G00_eV": 2.7,
        "cleave_gT_eV_per_K": 7.0e-5,
        "cleave_sigc0_GPa": 5.7,
        "emit_G00_eV": 3.1,
        "emit_gT_eV_per_K": 2.4e-3,
        "emit_sigc0_GPa": 4.4,
    }])


def run(dK: float):
    return atlas.evaluate_candidates(
        candidate(),
        shape_row(),
        [300.0, 700.0, 1200.0],
        Kdot=0.005,
        dK=dK,
        Kmax=50.0,
        nu0_c=1.0e12,
        nu0_e=1.0e11,
        m_hits=3.0,
        tau_c=1.0e-6,
        Tref_K=481.33,
        floor_min_eV=1.0e-4,
        floor_max_frac=0.95,
        progress_every=0,
    )


def test_exp_floor_decreases_with_stress():
    T = np.array([[700.0]])
    G00 = np.array([[2.0]])
    gT = np.array([[0.0]])
    sigc = np.array([[4.0]])
    low = atlas.exp_floor_barrier_eV(
        0.0,
        T,
        G00,
        gT,
        sigc,
        0.0,
        1.0,
        1.2,
        0.02,
        481.33,
        1.0e-4,
        0.95,
    )
    high = atlas.exp_floor_barrier_eV(
        10.0e9,
        T,
        G00,
        gT,
        sigc,
        0.0,
        1.0,
        1.2,
        0.02,
        481.33,
        1.0e-4,
        0.95,
    )
    assert float(high[0, 0]) < float(low[0, 0])


def test_analytic_first_passage_is_finite_and_emission_is_bounded():
    result = run(0.05)
    assert np.all(np.isfinite(result["Kc"]))
    assert np.all(result["Kc"] > 0.0)
    assert np.all(result["source_fraction_at_Kc"] >= 0.0)
    assert np.all(result["source_fraction_at_Kc"] <= 1.0)


def test_K_increment_refinement_is_stable():
    coarse = run(0.10)["Kc"]
    fine = run(0.02)["Kc"]
    assert np.nanmax(np.abs(coarse - fine)) < 0.15
