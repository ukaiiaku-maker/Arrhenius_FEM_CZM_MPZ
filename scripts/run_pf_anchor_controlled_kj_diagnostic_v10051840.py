"""Reproducible controlled-K_J diagnostic against a frozen PF anchor case.

This turns the first real controlled-K_J diagnostic run this workspace
performed (Peak, 1000 K, theta=0, rate1x, seed 8666 -- see
CLAUDE_PROGRESS.md's URGENT HANDOFF, item 9, and FEM_PF_PARITY_SCORECARD.md)
into a reusable, provenance-recording script instead of an ad hoc
`python3 -c "..."` invocation.

Scope -- read before using this script's output as evidence
--------------------------------------------------------------------------
This is a **Gate 1 (Numerical correctness) diagnostic only**. It calls
`sharp_front.main()` directly, bypassing the kernel-gated production
entry point: there is no signed shielding kernel, no persistent-site MPZ
engine, and no full-field bulk Peierls-Taylor closure in this
configuration (`front_state_model=legacy_scalar`). It proves the
safeguarded J/K_J-target controller can track a *prescribed* real PF
K_J(t) trajectory with the real solver, production-scale mesh, and the
correct Peak-class barrier row. It does **not** by itself demonstrate
FEM/PF physical correspondence -- see FEM_PF_PARITY_SCORECARD.md's
"Status summary" for why, and for what would be required to call it more
than that (extending it to first passage, or a native-loading run with
no PF target at all).

Provenance and fail-closed behavior
--------------------------------------------------------------------------
All PF inputs are resolved through
`arrhenius_fracture.pf_theta0_frozen_reference_v10051840.resolve_frozen_pf_reference`,
which never falls back to the live, externally-mutable PF repository --
only a local, frozen `reference_inputs/...` directory verified against a
recorded SHA-256 contract. If that directory is missing either file, or
either file's content no longer matches the recorded hash, resolution
raises immediately (fail closed) -- this happens even in `--dry-run`
mode, since it is a statement about input provenance, not about whether
the FEM solver is invoked.

Every run directory receives a `provenance_manifest.json` (git commit,
package version, resolved frozen-reference hashes, full assembled CLI
argv, anchor id, timestamp) written before the solver runs, so a run's
inputs remain auditable independent of this script's source history.
"""
from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(REPO_ROOT))

from arrhenius_fracture.pf_theta0_frozen_reference_v10051840 import (  # noqa: E402
    FrozenPFReference,
    resolve_frozen_pf_reference,
)

PACKAGE_NAME = "arrhenius-fem-czm"

# Peak-class cleavage/emission barrier row (candidate v913_zeroD_sobol_0242980,
# option v913_paper_peak01_0242980_persistent_sites, PF_v10.2.22_exact_persistent_site
# registry). Recorded verbatim from PF_REFERENCE_REGENERATION_CONTRACT.md's
# "Selected-row barrier/kinetic parameters" block -- this is a material-row
# property, not specific to any one PF run/campaign, so the same row applies
# to every temperature/rate anchor drawn from the Peak class.
PEAK_BARRIER_ROW = {
    "cleave_Tref_K": 481.33,
    "cleave_G00_eV": 4.011803912930191,
    "cleave_gT_eV_per_K": 0.0081468212408944,
    "cleave_sigc0_GPa": 7.349638969637454,
    "cleave_sT_GPa_per_K": -0.0012099820682778,
    "cleave_exp_a": 1.3227705595549195,
    "cleave_exp_n": 1.2558047282509506,
    "cleave_floor_frac": 0.0184513317135146,
    "emit_Tref_K": 481.33,
    "emit_G00_eV": 2.0446315124630927,
    "emit_gT_eV_per_K": 0.000941995373927,
    "emit_sigc0_GPa": 7.205215461552143,
    "emit_sT_GPa_per_K": 0.0013325903080403,
    "emit_exp_a": 0.0858719444833695,
    "emit_exp_n": 1.1423492641188204,
    "emit_floor_frac": 0.0058420097608974,
}


@dataclass(frozen=True)
class AnchorSpec:
    """Everything needed to reproduce one controlled-K_J diagnostic run."""

    anchor_id: str
    reference_dir: str
    expected_steps_csv_sha256: str
    expected_kernel_sha256: str
    temperature_K: float
    theta_deg: float
    dt_s: float
    n_stagger: int
    nx: int
    ny: int
    tip_h_fine_m: float
    tip_ratio: float
    da_phys_m: float
    dU_m: float
    front_state_model: str
    barrier_row: dict = field(default_factory=lambda: dict(PEAK_BARRIER_ROW))
    controller_rel_tol: float = 1.0e-3
    controller_max_iterations: int = 25
    default_steps: int = 40
    print_every: int = 4


ANCHOR_REGISTRY: dict[str, AnchorSpec] = {
    "peak_1000K_rate1x_seed8666_v10230": AnchorSpec(
        anchor_id="peak_1000K_rate1x_seed8666_v10230",
        reference_dir=str(
            REPO_ROOT
            / "reference_inputs/pf_final_v10_2_30/rate1x/"
            "v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666"
        ),
        expected_steps_csv_sha256=(
            "50153f1b93dee23407ec7971658fc8da917629c14e4ff3290b7a7a0c88449faa"
        ),
        expected_kernel_sha256=(
            "d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3"
        ),
        temperature_K=1000.0,
        theta_deg=0.0,
        dt_s=8.4,
        n_stagger=2,
        nx=36,
        ny=72,
        tip_h_fine_m=1.0e-6,
        tip_ratio=1.20,
        da_phys_m=5.0e-6,
        dU_m=5.0e-8,
        front_state_model="legacy_scalar",
    ),
}


def get_anchor(anchor_id: str) -> AnchorSpec:
    try:
        return ANCHOR_REGISTRY[anchor_id]
    except KeyError as exc:
        raise KeyError(
            f"unknown anchor id {anchor_id!r}; known anchors: "
            f"{sorted(ANCHOR_REGISTRY)}. Add a new AnchorSpec to "
            "ANCHOR_REGISTRY rather than constructing CLI args ad hoc."
        ) from exc


def resolve_anchor_reference(
    anchor: AnchorSpec, *, reference_dir_override: str | Path | None = None
) -> FrozenPFReference:
    """Fail-closed resolution of this anchor's frozen PF inputs.

    Raises FileNotFoundError/ValueError if the frozen directory is
    missing either file or either file's SHA-256 no longer matches the
    recorded contract -- never falls back to a live PF repository path.
    """
    return resolve_frozen_pf_reference(
        reference_dir_override if reference_dir_override is not None else anchor.reference_dir,
        expected_steps_csv_sha256=anchor.expected_steps_csv_sha256,
        expected_kernel_sha256=anchor.expected_kernel_sha256,
    )


def build_argv(
    anchor: AnchorSpec,
    frozen_ref: FrozenPFReference,
    out_dir: str | Path,
    *,
    steps: int | None = None,
) -> list[str]:
    """Assemble the exact `sharp_front` CLI argv for this anchor's diagnostic.

    Every value here is explicit (even where it matches `sharp_front`'s
    own default) so a future change to that default cannot silently
    change what this frozen diagnostic reproduces.
    """
    argv: list[str] = [
        "--out", str(out_dir),
        "--temperatures", str(anchor.temperature_K),
        "--steps", str(steps if steps is not None else anchor.default_steps),
        "--dt", str(anchor.dt_s),
        "--dU", str(anchor.dU_m),
        "--n-stagger", str(anchor.n_stagger),
        "--nx", str(anchor.nx),
        "--ny", str(anchor.ny),
        "--tip-h-fine", str(anchor.tip_h_fine_m),
        "--tip-ratio", str(anchor.tip_ratio),
        "--da-phys", str(anchor.da_phys_m),
        "--crystal-theta-deg", str(anchor.theta_deg),
        "--front-state-model", anchor.front_state_model,
        "--print-every", str(anchor.print_every),
        "--no-plots",
    ]
    for key, val in anchor.barrier_row.items():
        flag = "--" + key.replace("_", "-")
        argv.extend([flag, str(val)])
    argv.extend([
        "--pf-kj-target-csv", frozen_ref.steps_csv_path,
        "--pf-kj-target-rel-tol", str(anchor.controller_rel_tol),
        "--pf-kj-target-max-iterations", str(anchor.controller_max_iterations),
    ])
    return argv


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def _git_dirty() -> bool:
    try:
        out = subprocess.check_output(
            ["git", "status", "--short"], cwd=REPO_ROOT, text=True
        )
        return bool(out.strip())
    except Exception:
        return True


def _package_version() -> str:
    try:
        return version(PACKAGE_NAME)
    except PackageNotFoundError:
        return "unknown"


def build_provenance_manifest(
    anchor: AnchorSpec, frozen_ref: FrozenPFReference, argv: list[str]
) -> dict:
    return {
        "schema": "pf_anchor_controlled_kj_diagnostic_provenance_v10_0_5_18_4_0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "anchor_id": anchor.anchor_id,
        "git_commit": _git_commit(),
        "git_working_tree_dirty": _git_dirty(),
        "package_version": _package_version(),
        "frozen_reference": {
            "reference_dir": frozen_ref.reference_dir,
            "steps_csv_path": frozen_ref.steps_csv_path,
            "steps_csv_sha256": frozen_ref.steps_csv_sha256,
            "kernel_path": frozen_ref.kernel_path,
            "kernel_sha256": frozen_ref.kernel_sha256,
        },
        "argv": argv,
        "scope_note": (
            "Gate 1 / Numerical correctness only -- no signed kernel, no "
            "persistent-site MPZ engine, no full-field bulk Peierls-Taylor "
            "closure. See FEM_PF_PARITY_SCORECARD.md before citing this run "
            "as physical-correspondence evidence."
        ),
    }


def _default_out_dir(anchor: AnchorSpec) -> Path:
    return REPO_ROOT / "runs" / "anchor_diagnostics" / f"{anchor.anchor_id}_script"


def run_diagnostic(
    anchor_id: str,
    *,
    out_dir: str | Path | None = None,
    steps: int | None = None,
    reference_dir_override: str | Path | None = None,
    dry_run: bool = False,
    force: bool = False,
) -> dict:
    """Resolve inputs, assemble argv, write provenance, and (unless
    `dry_run`) invoke the real FEM solver. Returns the provenance manifest.
    """
    anchor = get_anchor(anchor_id)
    frozen_ref = resolve_anchor_reference(anchor, reference_dir_override=reference_dir_override)

    resolved_out = Path(out_dir) if out_dir is not None else _default_out_dir(anchor)
    if resolved_out.exists() and any(resolved_out.iterdir()) and not force:
        raise FileExistsError(
            f"output directory {resolved_out} already exists and is non-empty. "
            "Refusing to overwrite a prior run's evidence -- pass a different "
            "--out, or --force if you specifically intend to replace it."
        )
    resolved_out.mkdir(parents=True, exist_ok=True)

    argv = build_argv(anchor, frozen_ref, resolved_out, steps=steps)
    manifest = build_provenance_manifest(anchor, frozen_ref, argv)

    with open(resolved_out / "provenance_manifest.json", "w") as fp:
        json.dump(manifest, fp, indent=2, sort_keys=True)

    if dry_run:
        return manifest

    from arrhenius_fracture import sharp_front

    sharp_front.main(argv)
    return manifest


def _build_cli_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--anchor",
        default="peak_1000K_rate1x_seed8666_v10230",
        choices=sorted(ANCHOR_REGISTRY),
        help="Registered anchor case id.",
    )
    p.add_argument("--out", default=None, help="Output directory (default: derived from anchor id).")
    p.add_argument("--steps", type=int, default=None, help="Override the anchor's default step count.")
    p.add_argument(
        "--reference-dir",
        default=None,
        help="Override the frozen reference directory (advanced/testing use only).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve inputs, assemble argv, and write the provenance manifest without running the FEM solver.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Allow writing into a non-empty output directory (overwrites are otherwise refused).",
    )
    return p


def main(argv: list[str] | None = None) -> dict:
    args = _build_cli_parser().parse_args(argv)
    manifest = run_diagnostic(
        args.anchor,
        out_dir=args.out,
        steps=args.steps,
        reference_dir_override=args.reference_dir,
        dry_run=args.dry_run,
        force=args.force,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


if __name__ == "__main__":
    main()
