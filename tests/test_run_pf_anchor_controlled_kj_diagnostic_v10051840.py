from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest


def _module():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "run_pf_anchor_controlled_kj_diagnostic_v10051840.py"
    )
    spec = importlib.util.spec_from_file_location(
        "run_pf_anchor_controlled_kj_diagnostic_v10051840", path
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    # Register before exec: the module's dataclasses use `from __future__
    # import annotations` (postponed evaluation), and dataclass's own type
    # resolution looks the defining module up via sys.modules[cls.__module__]
    # -- it must already be registered there when exec_module runs the class
    # bodies, or that lookup returns None and dataclass() raises.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_frozen_reference(tmp_path: Path, csv_bytes: bytes, kernel_bytes: bytes) -> Path:
    ref_dir = tmp_path / "frozen_ref"
    ref_dir.mkdir()
    (ref_dir / "steps_1000K.csv").write_bytes(csv_bytes)
    (ref_dir / "family.json").write_bytes(kernel_bytes)
    return ref_dir


def test_known_anchor_is_registered():
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    assert anchor.temperature_K == 1000.0
    assert anchor.theta_deg == 0.0
    assert anchor.nx == 36 and anchor.ny == 72
    assert anchor.front_state_model == "legacy_scalar"
    # Barrier row must be the exact Peak-class row from
    # PF_REFERENCE_REGENERATION_CONTRACT.md, not a placeholder/default.
    assert anchor.barrier_row["cleave_G00_eV"] == pytest.approx(4.011803912930191)
    assert anchor.barrier_row["emit_G00_eV"] == pytest.approx(2.0446315124630927)


def test_unknown_anchor_raises_with_known_ids_listed():
    m = _module()
    with pytest.raises(KeyError, match="unknown anchor id"):
        m.get_anchor("does_not_exist")


def test_build_argv_contains_expected_flags_and_values(tmp_path: Path):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    frozen_ref = m.FrozenPFReference(
        reference_dir=str(tmp_path),
        steps_csv_path=str(tmp_path / "steps_1000K.csv"),
        steps_csv_sha256="deadbeef",
        kernel_path=str(tmp_path / "family.json"),
        kernel_sha256="feedface",
    )
    argv = m.build_argv(anchor, frozen_ref, tmp_path / "out", steps=7)

    def value_after(flag):
        return argv[argv.index(flag) + 1]

    assert value_after("--steps") == "7"
    assert value_after("--nx") == "36"
    assert value_after("--ny") == "72"
    assert value_after("--dt") == "8.4"
    assert value_after("--tip-h-fine") == "1e-06"
    assert value_after("--tip-ratio") == "1.2"
    assert value_after("--da-phys") == "5e-06"
    assert value_after("--crystal-theta-deg") == "0.0"
    assert value_after("--front-state-model") == "legacy_scalar"
    assert value_after("--cleave-G00-eV") == str(anchor.barrier_row["cleave_G00_eV"])
    assert value_after("--emit-G00-eV") == str(anchor.barrier_row["emit_G00_eV"])
    assert value_after("--pf-kj-target-csv") == frozen_ref.steps_csv_path
    assert value_after("--pf-kj-target-rel-tol") == str(anchor.controller_rel_tol)
    assert value_after("--pf-kj-target-max-iterations") == str(anchor.controller_max_iterations)
    assert "--no-plots" in argv


def test_build_argv_defaults_steps_from_anchor_when_not_overridden(tmp_path: Path):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    frozen_ref = m.FrozenPFReference(
        reference_dir=str(tmp_path),
        steps_csv_path=str(tmp_path / "steps_1000K.csv"),
        steps_csv_sha256="deadbeef",
        kernel_path=str(tmp_path / "family.json"),
        kernel_sha256="feedface",
    )
    argv = m.build_argv(anchor, frozen_ref, tmp_path / "out")
    assert argv[argv.index("--steps") + 1] == str(anchor.default_steps)


def test_argv_is_accepted_by_the_real_sharp_front_parser(tmp_path: Path):
    """The assembled argv must be valid input to sharp_front's own argparse
    (not just internally self-consistent) -- this is the actual CLI surface
    the diagnostic depends on."""
    from arrhenius_fracture import sharp_front

    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    frozen_ref = m.FrozenPFReference(
        reference_dir=str(tmp_path),
        steps_csv_path=str(tmp_path / "steps_1000K.csv"),
        steps_csv_sha256="deadbeef",
        kernel_path=str(tmp_path / "family.json"),
        kernel_sha256="feedface",
    )
    (tmp_path / "steps_1000K.csv").write_text("placeholder")
    argv = m.build_argv(anchor, frozen_ref, tmp_path / "out", steps=2)
    args = sharp_front._build_parser().parse_args(argv)
    assert args.nx == 36 and args.ny == 72
    assert args.cleave_G00_eV == pytest.approx(anchor.barrier_row["cleave_G00_eV"])
    assert args.pf_kj_target_csv == frozen_ref.steps_csv_path
    assert args.front_state_model == "legacy_scalar"


def test_resolve_anchor_reference_fails_closed_on_hash_mismatch(tmp_path: Path):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    ref_dir = _write_frozen_reference(tmp_path, b"tampered csv", b"tampered kernel")
    with pytest.raises(ValueError, match="Refusing to proceed"):
        m.resolve_anchor_reference(anchor, reference_dir_override=ref_dir)


def test_resolve_anchor_reference_fails_closed_on_missing_files(tmp_path: Path):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        m.resolve_anchor_reference(anchor, reference_dir_override=empty_dir)


def test_resolve_anchor_reference_succeeds_when_hashes_match(tmp_path: Path):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    csv_bytes = hashlib.sha256(b"").hexdigest().encode()  # arbitrary content
    kernel_bytes = b"kernel content"
    ref_dir = tmp_path / "matching_ref"
    ref_dir.mkdir()
    (ref_dir / "steps_1000K.csv").write_bytes(csv_bytes)
    (ref_dir / "family.json").write_bytes(kernel_bytes)

    # Build an anchor whose expected hashes match this tmp content, without
    # touching the module-level default anchor's real recorded hashes.
    import dataclasses

    patched = dataclasses.replace(
        anchor,
        expected_steps_csv_sha256=hashlib.sha256(csv_bytes).hexdigest(),
        expected_kernel_sha256=hashlib.sha256(kernel_bytes).hexdigest(),
    )
    resolved = m.resolve_anchor_reference(patched, reference_dir_override=ref_dir)
    assert resolved.steps_csv_sha256 == hashlib.sha256(csv_bytes).hexdigest()
    assert resolved.kernel_sha256 == hashlib.sha256(kernel_bytes).hexdigest()


def test_provenance_manifest_contains_required_fields(tmp_path: Path):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    frozen_ref = m.FrozenPFReference(
        reference_dir=str(tmp_path),
        steps_csv_path=str(tmp_path / "steps_1000K.csv"),
        steps_csv_sha256="deadbeef",
        kernel_path=str(tmp_path / "family.json"),
        kernel_sha256="feedface",
    )
    argv = m.build_argv(anchor, frozen_ref, tmp_path / "out", steps=3)
    manifest = m.build_provenance_manifest(anchor, frozen_ref, argv)

    assert manifest["anchor_id"] == anchor.anchor_id
    assert manifest["frozen_reference"]["steps_csv_sha256"] == "deadbeef"
    assert manifest["frozen_reference"]["kernel_sha256"] == "feedface"
    assert manifest["argv"] == argv
    assert "git_commit" in manifest
    assert "package_version" in manifest
    assert "generated_at_utc" in manifest
    assert "Gate 1" in manifest["scope_note"]


def test_dry_run_writes_manifest_and_does_not_invoke_the_solver(tmp_path: Path, monkeypatch):
    """--dry-run (via run_diagnostic(dry_run=True)) must resolve inputs and
    write provenance without ever importing/calling sharp_front.main --
    this is what lets this script's tests avoid launching a real FEM run."""
    m = _module()

    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    csv_bytes = b"synthetic steps csv"
    kernel_bytes = b"synthetic kernel"
    ref_dir = _write_frozen_reference(tmp_path, csv_bytes, kernel_bytes)

    import dataclasses

    patched = dataclasses.replace(
        anchor,
        anchor_id="test_anchor",
        expected_steps_csv_sha256=hashlib.sha256(csv_bytes).hexdigest(),
        expected_kernel_sha256=hashlib.sha256(kernel_bytes).hexdigest(),
    )
    monkeypatch.setitem(m.ANCHOR_REGISTRY, "test_anchor", patched)

    called = {"invoked": False}

    class _ExplodingSharpFront:
        @staticmethod
        def main(argv):
            called["invoked"] = True
            raise AssertionError("sharp_front.main must not run in dry-run mode")

    monkeypatch.setitem(
        __import__("sys").modules, "arrhenius_fracture.sharp_front", _ExplodingSharpFront
    )

    out_dir = tmp_path / "out"
    manifest = m.run_diagnostic(
        "test_anchor",
        out_dir=out_dir,
        reference_dir_override=ref_dir,
        dry_run=True,
    )

    assert called["invoked"] is False
    assert (out_dir / "provenance_manifest.json").is_file()
    on_disk = json.loads((out_dir / "provenance_manifest.json").read_text())
    assert on_disk["anchor_id"] == "test_anchor"
    assert manifest["anchor_id"] == "test_anchor"


def test_run_diagnostic_fails_closed_before_touching_output_dir_on_bad_reference(
    tmp_path: Path,
):
    """A tampered/missing frozen reference must be caught before any output
    directory or provenance manifest is created -- provenance must never be
    written for a run whose inputs could not be verified."""
    m = _module()
    ref_dir = _write_frozen_reference(tmp_path, b"tampered", b"tampered")
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Refusing to proceed"):
        m.run_diagnostic(
            "peak_1000K_rate1x_seed8666_v10230",
            out_dir=out_dir,
            reference_dir_override=ref_dir,
            dry_run=True,
        )
    assert not out_dir.exists()


def test_run_diagnostic_refuses_to_overwrite_nonempty_output_dir(tmp_path: Path, monkeypatch):
    m = _module()
    anchor = m.get_anchor("peak_1000K_rate1x_seed8666_v10230")
    csv_bytes = b"synthetic steps csv"
    kernel_bytes = b"synthetic kernel"
    ref_dir = _write_frozen_reference(tmp_path, csv_bytes, kernel_bytes)

    import dataclasses

    patched = dataclasses.replace(
        anchor,
        expected_steps_csv_sha256=hashlib.sha256(csv_bytes).hexdigest(),
        expected_kernel_sha256=hashlib.sha256(kernel_bytes).hexdigest(),
    )
    monkeypatch.setitem(m.ANCHOR_REGISTRY, "test_anchor2", patched)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "preexisting_evidence.json").write_text("{}")

    with pytest.raises(FileExistsError, match="non-empty"):
        m.run_diagnostic(
            "test_anchor2",
            out_dir=out_dir,
            reference_dir_override=ref_dir,
            dry_run=True,
        )


def test_cli_dry_run_end_to_end_against_real_frozen_anchor(tmp_path: Path):
    """End-to-end smoke of the actual CLI entry point against the real,
    committed frozen anchor bundle -- still --dry-run, so this never
    launches the FEM solver, only exercises resolution + argv assembly +
    provenance writing exactly as a real invocation would."""
    m = _module()
    out_dir = tmp_path / "cli_out"
    manifest = m.main(["--dry-run", "--out", str(out_dir)])
    assert manifest["anchor_id"] == "peak_1000K_rate1x_seed8666_v10230"
    assert (out_dir / "provenance_manifest.json").is_file()
