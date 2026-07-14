from pathlib import Path

import pandas as pd

from prepare_mpz_v9_3_pt_input import prepare_input


def test_prepare_input_restores_fixed_exp_floor_shape_columns(tmp_path: Path):
    atlas = pd.DataFrame([
        {
            "candidate_id": "weakT_r0p005_0000001",
            "region": "weakT_intrinsic",
            "shape_family": "weakT",
            "Kdot_MPa_sqrt_m_per_s": 0.005,
            "emit_G00_eV": 1.1,
            "emit_gT_eV_per_K": 0.002,
            "emit_sigc0_GPa": 2.0,
            "refined_Kc_T300": 15.0,
            "refined_Kc_T1200": 15.1,
        }
    ])
    materials = pd.DataFrame([
        {
            "analytic_candidate_id": "weakT_r0p005_0000001",
            "analytic_region": "weakT_intrinsic",
            "analytic_shape_family": "weakT",
            "emit_G00_eV": 99.0,
            "emit_gT_eV_per_K": 99.0,
            "emit_sigc0_GPa": 99.0,
            "emit_exp_a": 0.5,
            "emit_exp_n": 0.85,
            "emit_floor_frac": 0.03,
            "emit_sT_GPa_per_K": 0.0,
            "emit_Tref_K": 481.33,
        }
    ])
    ap = tmp_path / "atlas.csv"
    mp = tmp_path / "materials.csv"
    op = tmp_path / "joined.csv"
    atlas.to_csv(ap, index=False)
    materials.to_csv(mp, index=False)

    joined = prepare_input(ap, mp, op)

    assert op.exists()
    assert joined.loc[0, "emit_exp_a"] == 0.5
    assert joined.loc[0, "emit_exp_n"] == 0.85
    assert joined.loc[0, "emit_floor_frac"] == 0.03
    # Sampled intrinsic coordinates from the refined atlas remain authoritative.
    assert joined.loc[0, "emit_G00_eV"] == 1.1
    assert joined.loc[0, "emit_gT_eV_per_K"] == 0.002
    assert joined.loc[0, "emit_sigc0_GPa"] == 2.0
