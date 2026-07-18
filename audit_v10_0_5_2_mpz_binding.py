#!/usr/bin/env python3
"""Audit the authoritative v9.11 MPZ parser/factory binding for v10.0.5.2.

The v9.11 wrapper consumes ``--mpz-n-bins`` before the inner sharp-front parser
runs.  Consequently ``run_args.json`` contains the inner compatibility default
and is not evidence for the active moving-process-zone discretization.  This
audit replays the exact outer parser and factory configuration path and, for a
legacy completed output, ties that replay to the guarded runner source at the
recorded commit plus an explicit operator attestation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

from arrhenius_fracture import mixed_mode_first_passage_v9_11 as v911
from arrhenius_fracture.pf_equivalent_material_manifest import pf_manifest_path

PROVENANCE_NAME = "mpz_configuration_provenance_v10_0_5_2.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(path: Path) -> Any:
    return json.loads(path.read_text())


def _git_show(commit: str, path: str) -> str:
    result = subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def _historical_source_contract(commit: str, expected_bins: int) -> dict[str, Any]:
    runner = _git_show(commit, "run_v10_0_5_2_DBTT_700K_100um_gate.sh")
    v1003 = _git_show(commit, "arrhenius_fracture/mode_i_first_passage_v10_0_3_progressive.py")
    v911 = _git_show(commit, "arrhenius_fracture/mixed_mode_first_passage_v9_11.py")
    parameterization = _git_show(commit, "arrhenius_fracture/mpz_parameterization_v911.py")
    factory = _git_show(commit, "arrhenius_fracture/kinetic_campaign_czm_v10052.py")

    checks = {
        "guarded_runner_default": f"MPZ_N_BINS=${{MPZ_N_BINS:-{expected_bins}}}" in runner,
        "guarded_runner_rejects_other_bin_count": (
            'if [[ "$MPZ_N_BINS" != "200" ]]' in runner
            if expected_bins == 200 else True
        ),
        "guarded_runner_passes_outer_option": '--mpz-n-bins "$MPZ_N_BINS"' in runner,
        "v1003_preserves_explicit_outer_option": (
            '_option_value(remaining, "--mpz-n-bins") is None' in v1003
            and 'remaining.extend(["--mpz-n-bins", "200"])' in v1003
        ),
        "v911_outer_parser_owns_bin_option": (
            'p.add_argument("--mpz-n-bins", type=int, default=200)' in v911
        ),
        "v911_passes_outer_namespace_to_factory": (
            "_engine_factory(original_build, context, mm, row)" in v911
        ),
        "parameterization_assigns_outer_bin_count": (
            'cfg.n_bins = int(getattr(args, "mpz_n_bins", 200))' in parameterization
        ),
        "v10052_factory_uses_outer_namespace": (
            "cfg = v911.build_mpz_config(mm, row)" in factory
        ),
    }
    return {
        "commit": commit,
        "checks": checks,
        "verified": all(checks.values()),
    }


def replay_binding(material: str, length_um: float, n_bins: int) -> dict[str, Any]:
    manifest = pf_manifest_path(material)
    argv = [
        "--mpz-material-manifest", str(manifest),
        "--mpz-material-class", material,
        "--mpz-length-um", f"{float(length_um):.16g}",
        "--mpz-n-bins", str(int(n_bins)),
        "--target-traction-phase-deg", "0",
        "--reference-cleavage-shape", "1",
    ]
    mm, remaining = v911.parser().parse_known_args(argv)
    if remaining:
        raise RuntimeError(f"unexpected unparsed MPZ binding arguments: {remaining}")
    row = v911.load_selected_row(mm.mpz_material_manifest, mm.mpz_material_class)
    cfg = v911.build_mpz_config(mm, row)
    active_bins = int(cfg.n_bins)
    active_length_m = float(cfg.length_m)
    active_dx_m = active_length_m / active_bins
    verified = (
        active_bins == int(n_bins)
        and math.isclose(active_length_m, float(length_um) * 1.0e-6, rel_tol=0.0, abs_tol=1.0e-15)
    )
    return {
        "material_class": material,
        "parameter_manifest": str(manifest),
        "outer_parser_requested_mpz_n_bins": int(mm.mpz_n_bins),
        "outer_parser_requested_mpz_length_um": float(mm.mpz_length_um),
        "active_mpz_n_bins": active_bins,
        "active_mpz_length_m": active_length_m,
        "active_mpz_dx_m": active_dx_m,
        "active_mpz_source_bin_count": int(cfg.source_bin_count),
        "binding_replay_verified": verified,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--material", default="DBTT")
    parser.add_argument("--expected-length-um", type=float, default=100.0)
    parser.add_argument("--expected-mpz-bins", type=int, default=200)
    parser.add_argument("--run-commit", required=True)
    parser.add_argument("--guarded-runner-attested", action="store_true")
    args = parser.parse_args()

    root = args.root
    completion_path = root / "run_completion_v10_0_5_2.json"
    progressive_path = root / "kinetic_campaign_czm_progressive_2d_v10_0_3.json"
    channel_path = root / "parallel_channel_diagnostics_v10_0_5_2.json"
    run_args_path = root / "run_args.json"
    required = [completion_path, progressive_path, channel_path, run_args_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing completed output for MPZ provenance: " + ", ".join(missing))

    completion = dict(_load(completion_path))
    if completion.get("status") != "complete" or completion.get("run_completed_without_exception") is not True:
        raise RuntimeError("completion manifest does not certify a completed v10.0.5.2 invocation")

    binding = replay_binding(args.material, args.expected_length_um, args.expected_mpz_bins)
    historical = _historical_source_contract(args.run_commit, args.expected_mpz_bins)
    inner_args = dict(_load(run_args_path))
    verified = bool(
        args.guarded_runner_attested
        and binding["binding_replay_verified"]
        and historical["verified"]
    )

    payload = {
        "schema": "v10_0_5_2_mpz_configuration_provenance_v1",
        "point_release": "10.0.5.2",
        "provenance_method": "outer_parser_factory_binding_replay_with_guarded_runner_attestation",
        "runtime_state_directly_serialized_in_legacy_output": False,
        "guarded_runner_operator_attested": bool(args.guarded_runner_attested),
        "historical_source_contract": historical,
        **binding,
        "inner_sharp_front_run_args_mpz_n_bins_shadow": inner_args.get("mpz_n_bins"),
        "inner_shadow_is_authoritative": False,
        "inner_shadow_explanation": (
            "v9.11 consumes --mpz-n-bins in its outer parser; run_args.json is written "
            "from the inner sharp-front namespace and retains its compatibility default"
        ),
        "source_output_sha256": {
            path.name: _sha256(path)
            for path in (completion_path, progressive_path, channel_path, run_args_path)
        },
        "provenance_verified": verified,
    }
    path = root / PROVENANCE_NAME
    path.write_text(json.dumps(payload, indent=2, default=str))
    print(json.dumps(payload, indent=2, default=str))
    if not verified:
        raise RuntimeError("v10.0.5.2 MPZ binding provenance verification failed")
    print(f"V10.0.5.2 MPZ BINDING PROVENANCE VERIFIED: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
