"""v10.0.5.18: four-class FEM/CZM with stochastic cleavage and emission."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import mode_i_first_passage_v10_0_5_16_stochastic_pf_parity as _base
from .four_class_parameter_bridge_v100518 import (
    EXPECTED_OPTIONS,
    load_four_class_parameter_option,
)
from .persistent_site_registry_v100514 import ROWS
from .persistent_site_stochastic_emission_v100518 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518,
)

POINT_RELEASE = "10.0.5.18"
MODEL_ID = "FEM_CZM_four_class_stochastic_cleavage_emission_v10_0_5_18"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18.json"
DEFAULT_PARAMETER_SOURCE_ROOT = (
    Path(__file__).resolve().parents[1]
    / "runtime_inputs"
    / "v10_0_5_18_four_class"
)


class ProductionStochasticEmissionFrontEngineV100518(
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518
):
    """Production adapter preserving the legacy lambda_e diagnostic contract."""

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(
            K_cleave,
            K_emit,
            T,
            dt,
            metadata=metadata,
        )
        last = dict(self.mpz_state.last_emission or {})
        aggregate = np.asarray(
            last.get(
                "aggregate_hazard_final_by_system_s",
                last.get("aggregate_hazard_initial_by_system_s", np.zeros(2)),
            ),
            dtype=float,
        ).reshape(-1)
        out["lambda_e"] = float(np.sum(aggregate))
        out["lambda_e_semantics"] = (
            "instantaneous_sum_stochastic_signed_channel_aggregate_hazards"
        )
        return out



def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("--parameter-source-root", type=Path, default=None)
    parser.add_argument(
        "--parameter-option",
        choices=tuple(EXPECTED_OPTIONS),
        required=True,
    )
    return parser


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _option_value(argv: list[str], name: str) -> str | None:
    if name not in argv:
        return None
    index = argv.index(name)
    if index + 1 >= len(argv):
        return None
    return argv[index + 1]


def _replace_option(argv: list[str], name: str, value: str) -> None:
    while name in argv:
        index = argv.index(name)
        del argv[index : min(index + 2, len(argv))]
    argv.extend([name, str(value)])


def _single_temperature(argv: list[str]) -> float | str | None:
    raw = _option_value(argv, "--temperatures")
    if raw is None:
        return None
    tokens = str(raw).split()
    if len(tokens) == 1:
        try:
            return float(tokens[0])
        except ValueError:
            return tokens[0]
    return raw


def _rewrite_output_metadata(
    out: Path | None,
    audit: dict[str, Any],
    solver_args: list[str],
) -> None:
    if out is None:
        return
    out.mkdir(parents=True, exist_ok=True)
    theta_raw = _option_value(solver_args, "--crystal-theta-deg")
    theta_deg = None if theta_raw is None else float(theta_raw)
    audit = dict(audit)
    audit.update(
        {
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "theta_deg": theta_deg,
            "theta_is_runtime_input": True,
            "stochastic_cleavage_active": True,
            "stochastic_crack_advance_length_active": True,
            "stochastic_emission_active": True,
            "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
            "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
            "continuous_fractional_MPZ_translation": True,
            "microstructure_advance_precedes_cohesive_checkpoint": True,
        }
    )
    (out / TRANSFER_MANIFEST).write_text(
        json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n"
    )

    old_manifest = out / _base.PRODUCTION_MANIFEST
    if old_manifest.is_file():
        payload = json.loads(old_manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_18"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        payload["four_class_parameter_transfer"] = audit
        payload["parameter_option"] = audit["parameter_option"]
        payload["candidate_id"] = audit["candidate_id"]
        payload["material_class"] = audit["material_class"]
        payload["theta_deg"] = theta_deg
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(
            {
                "four_class_parameter_transfer_only": True,
                "parameter_values_manually_reconstructed": False,
                "parameter_values_refit": False,
                "PF_runtime_dependency": False,
                "theta_is_runtime_input": True,
                "stochastic_cleavage": True,
                "stochastic_crack_advance_length": True,
                "stochastic_emission": True,
                "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
                "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
                "emission_noise_added_to_K_J_stress_or_barriers": False,
                "continuous_fractional_MPZ_translation": True,
                "microstructure_advance_precedes_cohesive_checkpoint": True,
                "dynamic_tip_radius": True,
                "dynamic_front_width": True,
                "persistent_sites": True,
                "finite_source_inventory": False,
                "source_refresh": False,
                "explicit_recovery": False,
            }
        )
        payload["physics_contract"] = physics
        payload["front_engine"] = (
            ProductionStochasticEmissionFrontEngineV100518.audit_payload()
        )
        (out / PRODUCTION_MANIFEST).write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )
        old_manifest.unlink()

    old_selection = out / _base.SELECTION_MANIFEST
    if old_selection.is_file():
        selection = json.loads(old_selection.read_text())
        selection["schema"] = MODEL_ID
        selection["point_release"] = POINT_RELEASE
        selection["four_class_parameter_transfer"] = audit
        selection["parameter_option"] = audit["parameter_option"]
        selection["candidate_id"] = audit["candidate_id"]
        selection["material_class"] = audit["material_class"]
        selection["theta_deg"] = theta_deg
        policy = dict(selection.get("policy", {}) or {})
        policy.update(
            {
                "theta_is_runtime_input": True,
                "stochastic_cleavage": True,
                "stochastic_event_length": True,
                "stochastic_emission": True,
                "continuous_fractional_MPZ_translation": True,
                "PF_runtime_dependency": False,
            }
        )
        selection["policy"] = policy
        (out / SELECTION_MANIFEST).write_text(
            json.dumps(selection, indent=2, sort_keys=True, default=str) + "\n"
        )
        old_selection.unlink()

    engine_audit = ProductionStochasticEmissionFrontEngineV100518.audit_payload()
    cleavage = dict(engine_audit.get("stochastic_hazard", {}) or {})
    emission = dict(engine_audit.get("stochastic_emission", {}) or {})
    seed_payload = {
        "schema": "v10.0.5.18_domain_separated_stochastic_seed_manifest",
        "parameter_option": audit["parameter_option"],
        "candidate_id": audit["candidate_id"],
        "material_class": audit["material_class"],
        "temperature_K": _single_temperature(solver_args),
        "theta_deg": theta_deg,
        "cleavage_case_seed": cleavage.get("seed"),
        "emission_case_seed": emission.get("case_seed"),
        "emission_stream_seeds": emission.get("stream_seeds"),
        "emission_stream_semantics": [
            "signed_emission_channel_0",
            "signed_emission_channel_1",
        ],
        "rng": "numpy_default_rng",
        "parallel_job_order_independent": True,
    }
    (out / SEED_MANIFEST).write_text(
        json.dumps(seed_payload, indent=2, sort_keys=True, default=str) + "\n"
    )


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    wrapper, remaining = _parser().parse_known_args(user_args)
    source_root = (
        DEFAULT_PARAMETER_SOURCE_ROOT
        if wrapper.parameter_source_root is None
        else wrapper.parameter_source_root
    ).expanduser().resolve()
    candidate, audit = load_four_class_parameter_option(
        source_root,
        wrapper.parameter_option,
    )
    out = _out_path(remaining)
    _replace_option(remaining, "--persistent-site-option", candidate.option_key)

    existed = candidate.option_key in ROWS
    previous = ROWS.get(candidate.option_key)
    if existed and previous != candidate:
        raise RuntimeError(
            f"runtime registry already contains a different row for {candidate.option_key}"
        )
    ROWS[candidate.option_key] = candidate
    saved_engine = _base.PersistentSitePFStochasticMovingTipFrontEngineV100516
    _base.PersistentSitePFStochasticMovingTipFrontEngineV100516 = (
        ProductionStochasticEmissionFrontEngineV100518
    )
    try:
        return _base.main(remaining)
    finally:
        _rewrite_output_metadata(out, audit, remaining)
        _base.PersistentSitePFStochasticMovingTipFrontEngineV100516 = saved_engine
        if existed:
            ROWS[candidate.option_key] = previous
        else:
            ROWS.pop(candidate.option_key, None)


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_PARAMETER_SOURCE_ROOT",
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "ProductionStochasticEmissionFrontEngineV100518",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
