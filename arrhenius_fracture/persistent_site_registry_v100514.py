"""Exact PF v10.2.22 persistent-site parameter registry for FEM/CZM parity."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from typing import Any

POINT_RELEASE = "10.0.5.14"
PARAMETER_SOURCE = "PF_v10.2.22_exact_persistent_site_registry"


@dataclass(frozen=True)
class PersistentSiteRowV100514:
    option_key: str
    candidate_id: str
    role: str
    emit_G00_eV: float
    emit_gT_eV_per_K: float
    peierls_H0_eV: float
    peierls_activation_entropy_kB: float
    taylor_H0_eV: float
    taylor_activation_entropy_kB: float
    taylor_corr_rho_c_m2: float
    rho_source0_m2: float
    source_refresh_length_um_provenance: float

    # Shared PF v10.2.22 contract.
    Tref_K: float = 481.33
    n_slip_channels: int = 2
    rho_forest_floor_m2: float = 5.0e12
    peierls_stress_fraction: float = 0.5773502691896258
    taylor_stress_fraction: float = 0.5773502691896258
    mobile_shield_fraction: float = 0.0
    source_recovery_rate_s: float = 0.0
    L_pz_um_recommended: float = 50.0
    n_bins_recommended: int = 80
    cleave_G00_eV: float = 3.412309742899366
    cleave_gT_eV_per_K: float = 0.0085569133975694
    cleave_sigc0_GPa: float = 4.178054527802145
    cleave_sT_GPa_per_K: float = 0.004399175167258
    cleave_exp_a: float = 0.4846029535103278
    cleave_exp_n: float = 2.115325382808986
    cleave_floor_frac: float = 0.0543870261570049
    emit_sigc0_GPa: float = 5.224677563435593
    emit_sT_GPa_per_K: float = -0.0045302835171498
    emit_exp_a: float = 0.2784130973036495
    emit_exp_n: float = 0.8357097011358764
    emit_floor_frac: float = 0.068706727639919
    peierls_exp_a: float = 1.14391738514239
    peierls_exp_n: float = 1.7640893795753347
    peierls_nu0_s: float = 1.0e12
    taylor_exp_a: float = 0.1601883464066985
    taylor_exp_n: float = 1.5292178333131945
    taylor_nu0_s: float = 1.0e11
    taylor_corr_scale: float = 1.3270040200432451
    source_sites_per_system_provenance: float = 141.0590567476921
    encounter_efficiency: float = 9.160246308716648
    retained_recovery_rate_s: float = 0.0
    c_blunt: float = 1.411283192139077
    recovery_nu0_s: float = 0.0
    reference_source_area_um2: float = 25.0
    reference_front_width_um: float = 10.0
    source_zone_length_um: float = 2.0
    minimum_front_width_um: float = 0.0

    def validate(self) -> "PersistentSiteRowV100514":
        if self.rho_source0_m2 <= 0.0:
            raise ValueError("rho_source0_m2 must be positive")
        if self.n_slip_channels != 2:
            raise ValueError("v10.0.5.14 requires exactly two reduced slip channels")
        for name in (
            "source_recovery_rate_s",
            "retained_recovery_rate_s",
            "recovery_nu0_s",
        ):
            if abs(float(getattr(self, name))) > 1.0e-30:
                raise ValueError(f"{name} must be zero for PF v10.2.22 parity")
        return self

    def barrier_row(self) -> dict[str, Any]:
        return {
            "Tref_K": self.Tref_K,
            "cleave_G00_eV": self.cleave_G00_eV,
            "cleave_gT_eV_per_K": self.cleave_gT_eV_per_K,
            "cleave_sigc0_GPa": self.cleave_sigc0_GPa,
            "cleave_sT_GPa_per_K": self.cleave_sT_GPa_per_K,
            "cleave_exp_a": self.cleave_exp_a,
            "cleave_exp_n": self.cleave_exp_n,
            "cleave_floor_frac": self.cleave_floor_frac,
            "emit_G00_eV": self.emit_G00_eV,
            "emit_gT_eV_per_K": self.emit_gT_eV_per_K,
            "emit_sigc0_GPa": self.emit_sigc0_GPa,
            "emit_sT_GPa_per_K": self.emit_sT_GPa_per_K,
            "emit_exp_a": self.emit_exp_a,
            "emit_exp_n": self.emit_exp_n,
            "emit_floor_frac": self.emit_floor_frac,
            "peierls_H0_eV": self.peierls_H0_eV,
            "peierls_activation_entropy_kB": self.peierls_activation_entropy_kB,
            "peierls_exp_a": self.peierls_exp_a,
            "peierls_exp_n": self.peierls_exp_n,
            "peierls_nu0_s": self.peierls_nu0_s,
            "taylor_H0_eV": self.taylor_H0_eV,
            "taylor_activation_entropy_kB": self.taylor_activation_entropy_kB,
            "taylor_exp_a": self.taylor_exp_a,
            "taylor_exp_n": self.taylor_exp_n,
            "taylor_nu0_s": self.taylor_nu0_s,
        }

    def fingerprint(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()


ROWS = {
    "v912_top1_peak_persistent_sites": PersistentSiteRowV100514(
        "v912_top1_peak_persistent_sites",
        "v912_targeted_local_peak_013476_0368",
        "Primary sharp peak candidate 0368",
        1.852910757251084,
        0.001915784010291,
        4.962303673103452,
        -8.27728333324194,
        0.5828757468677228,
        17.22111565526575,
        6.549380478082156e16,
        1.4115242646890916e16,
        5.126255443312889,
    ),
    "v912_peak_0314_persistent_sites": PersistentSiteRowV100514(
        "v912_peak_0314_persistent_sites",
        "v912_targeted_local_peak_013476_0314",
        "Broad-onset peak candidate 0314",
        1.595293372310698,
        0.0025668888360261,
        4.588103987649083,
        -15.637479491531847,
        0.4820600483573401,
        23.443809836171567,
        4.121272046109432e15,
        1.4854598439174714e16,
        43.39285168360589,
    ),
    "v912_peak_0162_persistent_sites": PersistentSiteRowV100514(
        "v912_peak_0162_persistent_sites",
        "v912_targeted_local_peak_013476_0162",
        "Lower-peak candidate 0162",
        1.7407055431976914,
        0.0025359431378543,
        4.918964824452996,
        -5.982995696365832,
        0.5551850934439608,
        23.317227921448648,
        1.2821927846870708e16,
        2.3391211664562664e16,
        9.174389554897278,
    ),
    "v912_peak_0118_persistent_sites": PersistentSiteRowV100514(
        "v912_peak_0118_persistent_sites",
        "v912_targeted_local_peak_005518_0118",
        "Conventional DBTT topology candidate 0118",
        4.654189503333349,
        -0.0020404944046833,
        0.1896856116048691,
        -12.12842506536634,
        0.6090008780360222,
        17.170955220237374,
        8.604296690405017e11,
        2.3012695321899512e16,
        13.825934647000269,
    ),
    "v912_plateau_0403_persistent_sites": PersistentSiteRowV100514(
        "v912_plateau_0403_persistent_sites",
        "v912_targeted_local_plateau_010759_0403",
        "Plateau control candidate 0403",
        1.605401717312634,
        0.0014598098032176,
        0.6528542017564176,
        7.838268833234905,
        0.3551186254143289,
        12.783966001842384,
        3.554616121142454e17,
        1.3345163074387614e16,
        30.332215564960894,
    ),
}


def select_persistent_site_row(option_key: str) -> PersistentSiteRowV100514:
    try:
        return ROWS[str(option_key)].validate()
    except KeyError as exc:
        raise KeyError(
            f"unknown persistent-site option {option_key!r}; allowed={sorted(ROWS)}"
        ) from exc


__all__ = [
    "POINT_RELEASE",
    "PARAMETER_SOURCE",
    "PersistentSiteRowV100514",
    "ROWS",
    "select_persistent_site_row",
]
