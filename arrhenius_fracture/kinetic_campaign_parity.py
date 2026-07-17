"""Isolated prescribed-K parity harness for ``kinetic_campaign_czm``.

This harness removes the FEM and cohesive feedback.  It is the Stage-A gate:
the CZM front state is driven by an explicit K/T/dt history and exports every
state quantity required for comparison with PF v10.1.7.1.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable

import numpy as np

from .config import EV_TO_J
from .kinetic_campaign_czm import (
    DevelopedStateDiagnosticCZMFrontEngine,
    KineticCampaignCZMConfig,
    apply_pf_manifest_to_mpz_config,
)
from .moving_process_zone import MovingProcessZoneConfig
from .pf_equivalent_material_manifest import (
    PF_SOURCE,
    load_material_manifest,
)


class _ManifestBarrierAdapter:
    barrier_kind = "exp_floor"

    def __init__(self, surface):
        self.surface = surface

    def G_barrier(self, stress_Pa, T_K, b):
        return self.surface.values_eV(stress_Pa, T_K) * EV_TO_J

    def diagnostics(self, stress_Pa, T_K, b):
        sigma = np.asarray(stress_Pa, dtype=float)
        G = self.surface.values_eV(sigma, T_K)
        eps = 1.0e5
        gp = self.surface.values_eV(sigma + eps, T_K)
        gm = self.surface.values_eV(np.maximum(sigma - eps, 0.0), T_K)
        derivative = (gp - gm) / (2.0 * eps) * 1.0e9
        return {
            "G_eV": G,
            "S_kB": np.zeros_like(G),
            "dG_dsigma_eV_per_GPa": derivative,
            "vstar_b3": np.zeros_like(G),
        }


def build_isolated_campaign_engine(
    material_class: str,
    *,
    parameter_source: str = PF_SOURCE,
    mpz_length_m: float = 100.0e-6,
    mpz_n_bins: int = 200,
    da_phys_m: float = 5.0e-6,
    G_Pa: float = 80.0e9,
    nu: float = 0.28,
    b_m: float = 2.48e-10,
    r0_m: float = 1.0e-6,
    kinetic_config: KineticCampaignCZMConfig | None = None,
) -> DevelopedStateDiagnosticCZMFrontEngine:
    manifest = load_material_manifest(
        material_class, parameter_source=parameter_source
    )
    kinetic = (kinetic_config or KineticCampaignCZMConfig()).validate()
    cfg = MovingProcessZoneConfig(
        length_m=float(mpz_length_m),
        n_bins=int(mpz_n_bins),
        n_systems=2,
    )
    apply_pf_manifest_to_mpz_config(cfg, manifest, kinetic)
    fcfg = SimpleNamespace(
        r0=float(r0_m),
        L_pz=float(mpz_length_m),
        c_blunt=float(manifest.c_blunt),
        nu0_c=float(manifest.cleavage.attempt_frequency_s),
        nu0_e=float(manifest.emission.attempt_frequency_s),
        m_hits=3.0,
        tau_c=1.0e-6,
        sigma_cap=0.0,
        dN_cap=math.inf,
        N_sat=math.inf,
        recover_k=0.0,
        k_shield=0.0,
        chi_shield=0.0,
        v_emb_b3=0.0,
        emb_sat_frac=1.0,
        beta_back=0.0,
        rho0=5.0e12,
        tau_B=0.0,
        da=float(da_phys_m),
        max_advances_per_step=1.0,
    )
    engine = DevelopedStateDiagnosticCZMFrontEngine(
        fcfg,
        _ManifestBarrierAdapter(manifest.cleavage),
        _ManifestBarrierAdapter(manifest.emission),
        float(G_Pa),
        float(nu),
        float(b_m),
        cfg,
        manifest,
        kinetic,
    )
    return engine


def _default_history() -> list[dict[str, float]]:
    rows = []
    for i, K in enumerate(np.linspace(2.0, 35.0, 200)):
        rows.append({
            "step": float(i),
            "K_open_MPa_sqrt_m": float(K),
            "K_cleave_MPa_sqrt_m": float(K),
            "temperature_K": 700.0,
            "dt_s": 1.0e-4,
            "weight_0": 1.0,
            "weight_1": 1.0,
        })
    return rows


def load_history(path: str | Path | None) -> list[dict[str, float]]:
    if path is None:
        return _default_history()
    with Path(path).open(newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    out = []
    for i, row in enumerate(rows):
        out.append({
            "step": float(row.get("step", i)),
            "K_open_MPa_sqrt_m": float(row["K_open_MPa_sqrt_m"]),
            "K_cleave_MPa_sqrt_m": float(
                row.get("K_cleave_MPa_sqrt_m", row["K_open_MPa_sqrt_m"])
            ),
            "temperature_K": float(row["temperature_K"]),
            "dt_s": float(row["dt_s"]),
            "weight_0": float(row.get("weight_0", 1.0)),
            "weight_1": float(row.get("weight_1", 1.0)),
        })
    return out


def run_history(engine, history: Iterable[dict[str, float]]) -> list[dict[str, Any]]:
    output = []
    for row in history:
        weights = np.array([row.get("weight_0", 1.0), row.get("weight_1", 1.0)])
        result = engine.integrate_kinetics(
            float(row["K_open_MPa_sqrt_m"]) * 1.0e6,
            float(row["K_cleave_MPa_sqrt_m"]) * 1.0e6,
            float(row["temperature_K"]),
            float(row["dt_s"]),
            system_weights=weights,
        )
        channels = engine.stress_channels(
            float(row["K_open_MPa_sqrt_m"]) * 1.0e6,
            float(row["K_cleave_MPa_sqrt_m"]) * 1.0e6,
            weights,
        )
        state = engine.mpz_state.diagnostics_campaign()
        output.append({
            **row,
            "B": float(engine.B),
            "micro_advance_total_m": float(engine.micro_advance_total_m),
            "checkpoint_advance_total_m": float(engine.checkpoint_advance_total_m),
            "n_advances": int(engine.n_adv),
            "time_s": float(engine.t),
            "fired": bool(result["fired"]),
            "dt_consumed_s": float(result["dt_consumed_s"]),
            "dt_unused_s": float(result["dt_unused_s"]),
            "internal_substeps": int(result["internal_substeps"]),
            "r_eff_m": float(engine.r_eff()),
            **channels,
            **state,
        })
        if result["fired"] and result["dt_unused_s"] > 0.0:
            break
    return output


def compare_records(
    reference: list[dict[str, Any]],
    candidate: list[dict[str, Any]],
    *,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-12,
) -> dict[str, Any]:
    keys = (
        "B",
        "micro_advance_total_m",
        "source_budget_remaining",
        "mobile_count",
        "retained_count",
        "sigma_emission_backstress_Pa",
        "K_shield_effective_Pa_sqrt_m",
        "r_eff_m",
        "cumulative_emitted",
        "cumulative_refreshed",
        "cumulative_trapped",
        "cumulative_released",
        "cumulative_recovered",
        "cumulative_escaped",
    )
    n = min(len(reference), len(candidate))
    comparisons = []
    passed = len(reference) == len(candidate)
    for key in keys:
        ref = np.array([float(reference[i][key]) for i in range(n)])
        cand = np.array([float(candidate[i][key]) for i in range(n)])
        ok = bool(np.allclose(ref, cand, rtol=rtol, atol=atol))
        passed = passed and ok
        comparisons.append({
            "field": key,
            "passed": ok,
            "max_abs_error": float(np.max(np.abs(ref - cand))) if n else None,
            "max_rel_error": float(
                np.max(np.abs(ref - cand) / np.maximum(np.abs(ref), atol))
            ) if n else None,
        })
    return {
        "schema": "pf_czm_isolated_front_parity_v10_0",
        "reference_count": len(reference),
        "candidate_count": len(candidate),
        "rtol": rtol,
        "atol": atol,
        "passed": passed,
        "fields": comparisons,
    }


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--material", choices=("ceramic", "weakT", "DBTT"), default="weakT")
    p.add_argument("--history-csv")
    p.add_argument("--out", required=True)
    p.add_argument("--reference-json")
    args = p.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    history = load_history(args.history_csv)
    engine = build_isolated_campaign_engine(args.material)
    records = run_history(engine, history)
    (out / "czm_isolated_front_records.json").write_text(
        json.dumps(records, indent=2, default=str)
    )
    with (out / "czm_isolated_front_records.csv").open("w", newline="") as handle:
        if records:
            writer = csv.DictWriter(handle, fieldnames=sorted(records[0]))
            writer.writeheader()
            writer.writerows(records)
    audit = {
        "schema": "kinetic_campaign_czm_isolated_stage_A_v10_0",
        "material": args.material,
        "model_audit": engine.audit_payload(),
        "record_count": len(records),
        "reference_comparison_performed": bool(args.reference_json),
    }
    if args.reference_json:
        reference = json.loads(Path(args.reference_json).read_text())
        comparison = compare_records(reference, records)
        audit["comparison"] = comparison
        if not comparison["passed"]:
            (out / "isolated_front_parity_audit.json").write_text(
                json.dumps(audit, indent=2, default=str)
            )
            raise SystemExit("PF/CZM isolated front parity failed")
    (out / "isolated_front_parity_audit.json").write_text(
        json.dumps(audit, indent=2, default=str)
    )
    print(
        f"ISOLATED KINETIC CZM {args.material}: records={len(records)} "
        f"advance_um={engine.checkpoint_advance_total_m * 1e6:.6f} "
        f"B={engine.B:.9f}"
    )
    return records


if __name__ == "__main__":
    main()
