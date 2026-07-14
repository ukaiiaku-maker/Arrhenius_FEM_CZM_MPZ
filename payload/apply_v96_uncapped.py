from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    if new in text:
        return
    if old not in text:
        raise SystemExit(f"patch anchor not found in {path}: {old[:140]!r}")
    p.write_text(text.replace(old, new, 1))


# Make scalar and array harmonic rates robust without artificial epsilon floors.
replace_once(
    "arrhenius_fracture/emission_derived_plasticity_v96.py",
    '''        out = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=float)
        aa, bb = np.broadcast_arrays(a, b)
        active = (aa > 0.0) & (bb > 0.0)
        out[active] = aa[active] * bb[active] / (aa[active] + bb[active])
        return out
''',
    '''        aa, bb = np.broadcast_arrays(a, b)
        active = (aa > 0.0) & (bb > 0.0)
        denom = aa + bb
        return np.where(active, aa * bb / np.where(active, denom, 1.0), 0.0)
''',
)

# Attach each exact historical reference to the nearest fully evaluated atlas
# first-passage curve. The match is a transparent proxy for K(T), not a claim
# that the barriers are identical.
replace_once(
    "search_mpz_v9_6_broad_dbtt_map.py",
    '''def normalize_canonical(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.copy()
    df["candidate_id"] = "canonical_" + df.target_class.astype(str)
    df["region"] = df.target_class.map({
        "ceramic": "ceramic_reference",
        "peak": "peak_reference",
        "weakT": "weakT_reference",
        "DBTT": "DBTT_reference",
    })
    df["candidate_source"] = "prior_first_passage_reference"
    return df
''',
    '''SURFACE_MATCH_COLUMNS = [
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac", "cleave_S_hs_kB",
]


def normalize_canonical(
    path: Path, intrinsic: pd.DataFrame, temperatures: list[float]
) -> pd.DataFrame:
    refs = pd.read_csv(path).copy()
    refs["candidate_id"] = "canonical_" + refs.target_class.astype(str)
    refs["region"] = refs.target_class.map({
        "ceramic": "ceramic_reference",
        "peak": "peak_reference",
        "weakT": "weakT_reference",
        "DBTT": "DBTT_reference",
    })
    refs["candidate_source"] = "prior_first_passage_reference"

    kc_names = []
    for T in temperatures:
        tag = f"{int(round(T))}"
        name = next(
            (c for c in (f"refined_Kc_T{tag}", f"Kc_T{tag}") if c in intrinsic),
            None,
        )
        if name is None:
            raise ValueError(f"atlas has no first-passage K column for {T:g} K")
        kc_names.append(name)
    pool = intrinsic.dropna(subset=kc_names).copy()
    common = [c for c in SURFACE_MATCH_COLUMNS if c in refs and c in pool]
    if not common or pool.empty:
        raise ValueError("cannot match canonical references to evaluated atlas rows")
    X = pool[common].to_numpy(float)
    scale = np.nanstd(X, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-12), scale, 1.0)
    for idx, ref in refs.iterrows():
        y = ref[common].to_numpy(float)
        distance = np.sqrt(np.nanmean(((X - y) / scale) ** 2, axis=1))
        j = int(np.nanargmin(distance))
        nearest = pool.iloc[j]
        for name in kc_names:
            refs.loc[idx, name] = float(nearest[name])
        refs.loc[idx, "canonical_kc_proxy_candidate_id"] = str(nearest.candidate_id)
        refs.loc[idx, "canonical_kc_proxy_distance"] = float(distance[j])
    return refs
''',
)
replace_once(
    "search_mpz_v9_6_broad_dbtt_map.py",
    '''    canonical = normalize_canonical(a.canonical)
    candidates = pd.concat([intrinsic, canonical], ignore_index=True, sort=False)
''',
    '''    temperatures = floats(a.temperatures)
    canonical = normalize_canonical(a.canonical, intrinsic, temperatures)
    candidates = pd.concat([intrinsic, canonical], ignore_index=True, sort=False)
''',
)
replace_once(
    "search_mpz_v9_6_broad_dbtt_map.py",
    '''    temperatures = floats(a.temperatures)
    targets = target_curve(a.targets, temperatures)
''',
    '''    targets = target_curve(a.targets, temperatures)
''',
)
replace_once(
    "search_mpz_v9_6_broad_dbtt_map.py",
    '''                    "region": row.get("region", ""),
                    **p.to_dict(),
                    "T_K": T,
''',
    '''                    "region": row.get("region", ""),
                    "canonical_kc_proxy_candidate_id": row.get("canonical_kc_proxy_candidate_id", ""),
                    "canonical_kc_proxy_distance": row.get("canonical_kc_proxy_distance", np.nan),
                    **p.to_dict(),
                    "T_K": T,
''',
)
replace_once(
    "search_mpz_v9_6_broad_dbtt_map.py",
    '''                    "region": row.get("region", ""),
                    **p.to_dict(),
                    "dbtt_proxy_score": score,
''',
    '''                    "region": row.get("region", ""),
                    "canonical_kc_proxy_candidate_id": row.get("canonical_kc_proxy_candidate_id", ""),
                    "canonical_kc_proxy_distance": row.get("canonical_kc_proxy_distance", np.nan),
                    **p.to_dict(),
                    "dbtt_proxy_score": score,
''',
)

# Activate v9.6 everywhere the production modules lazily import the PT model.
replace_once(
    "arrhenius_fracture/__init__.py",
    '''from .moving_process_zone_v95 import MovingProcessZoneState as _MPZStateV95
_mpz_v94.MovingProcessZoneState = _MPZStateV95
MovingProcessZoneState = _MPZStateV95

__version__ = '0.9.5'
''',
    '''from .moving_process_zone_v95 import MovingProcessZoneState as _MPZStateV95
_mpz_v94.MovingProcessZoneState = _MPZStateV95
MovingProcessZoneState = _MPZStateV95

# v9.6 removes the exploratory PT caps and saturation functions while retaining
# the v9.4 detailed-balance law and the v9.5 spatial MPZ state.
from .emission_derived_plasticity_v96 import (
    EmissionDerivedPeierlsTaylorModel as _PTModelV96,
)
_pt_v93.EmissionDerivedPeierlsTaylorModel = _PTModelV96

__version__ = '0.9.6'
''',
)
replace_once("pyproject.toml", 'version = "0.9.5"', 'version = "0.9.6"')
