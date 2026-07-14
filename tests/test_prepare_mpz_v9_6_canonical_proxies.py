from pathlib import Path

import pandas as pd

from prepare_mpz_v9_6_canonical_proxies import prepare


def test_prepare_attaches_complete_atlas_curve_and_records_match(tmp_path: Path):
    atlas = pd.DataFrame([
        {
            "candidate_id": "near",
            "emit_G00_eV": 2.2,
            "emit_gT_eV_per_K": 0.001,
            "emit_sigc0_GPa": 3.7,
            "emit_sT_GPa_per_K": 0.002,
            "emit_exp_a": 0.1,
            "emit_exp_n": 0.84,
            "emit_floor_frac": 0.035,
            "cleave_G00_eV": 2.8,
            "cleave_gT_eV_per_K": 0.0045,
            "cleave_sigc0_GPa": 4.1,
            "cleave_sT_GPa_per_K": 0.0001,
            "cleave_exp_a": 0.59,
            "cleave_exp_n": 1.2,
            "cleave_floor_frac": 0.0026,
            "cleave_S_hs_kB": 7.0,
            "refined_Kc_T300": 15.0,
            "refined_Kc_T700": 14.0,
            "refined_Kc_T900": 11.0,
            "refined_Kc_T1200": 8.0,
        },
        {
            "candidate_id": "far",
            "emit_G00_eV": 9.0,
            "emit_gT_eV_per_K": 0.02,
            "emit_sigc0_GPa": 9.0,
            "emit_sT_GPa_per_K": 0.01,
            "emit_exp_a": 1.5,
            "emit_exp_n": 2.0,
            "emit_floor_frac": 0.1,
            "cleave_G00_eV": 9.0,
            "cleave_gT_eV_per_K": 0.02,
            "cleave_sigc0_GPa": 9.0,
            "cleave_sT_GPa_per_K": 0.01,
            "cleave_exp_a": 1.5,
            "cleave_exp_n": 2.0,
            "cleave_floor_frac": 0.1,
            "cleave_S_hs_kB": -20.0,
            "refined_Kc_T300": 40.0,
            "refined_Kc_T700": 40.0,
            "refined_Kc_T900": 40.0,
            "refined_Kc_T1200": 40.0,
        },
    ])
    canonical = atlas.iloc[[0]].drop(
        columns=[c for c in atlas if c.startswith("refined_Kc_")]
    ).copy()
    canonical["target_class"] = "DBTT"
    canonical = canonical.drop(columns=["candidate_id"])

    ap = tmp_path / "atlas.csv"
    cp = tmp_path / "canonical.csv"
    op = tmp_path / "proxy.csv"
    atlas.to_csv(ap, index=False)
    canonical.to_csv(cp, index=False)

    result = prepare(ap, cp, [300, 700, 900, 1200], op)

    assert op.exists()
    assert result.loc[0, "canonical_kc_proxy_candidate_id"] == "near"
    assert result.loc[0, "refined_Kc_T300"] == 15.0
    assert result.loc[0, "refined_Kc_T1200"] == 8.0
    assert result.loc[0, "candidate_source"] == "prior_first_passage_reference"
