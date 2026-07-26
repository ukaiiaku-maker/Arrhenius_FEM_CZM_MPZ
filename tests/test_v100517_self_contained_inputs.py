from __future__ import annotations

from pathlib import Path

from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_17_paper_parameter_campaign as campaign,
)

RUNNER = Path("run_v10_0_5_17_paper7_300_1200K_200um_stochastic_campaign.sh")
VALIDATED = Path("run_v10_0_5_17_paper7_validated_campaign.sh")
VALIDATOR = Path("scripts/validate_v100517_paper_parameter_inputs.py")
VENDOR = Path("scripts/vendor_v100517_campaign_inputs.py")


def test_default_parameter_source_is_local_to_fem_checkout():
    expected = (
        Path(campaign.__file__).resolve().parents[1]
        / "runtime_inputs"
        / "v10_0_5_17_frozen_pf_inputs"
    )
    assert campaign.DEFAULT_PARAMETER_SOURCE_ROOT == expected


def test_production_launchers_have_no_live_pf_checkout_dependency():
    text = RUNNER.read_text() + "\n" + VALIDATED.read_text()
    for forbidden in (
        "PFROOT",
        "PF_BRANCH_REQUIRED",
        'git -C "$PFROOT"',
        "--pf-repo-root",
    ):
        assert forbidden not in text
    assert "PARAMETER_SOURCE_ROOT" in text
    assert "--parameter-source-root" in text
    assert "external_PF_runtime_dependency=false" in text


def test_original_physics_command_is_preserved():
    text = RUNNER.read_text()
    required = (
        "CLEAVAGE_HAZARD_MODE=exponential",
        "CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled",
        'CLEAVAGE_EVENT_MIN_FACTOR="$EVENT_MIN_FACTOR"',
        'CLEAVAGE_EVENT_MAX_FACTOR="$EVENT_MAX_FACTOR"',
        "--crack-backend adaptive_czm",
        "--crystal-aniso",
        "--crystal-compete",
        "--crystal-theta-deg 45",
        "--max-fronts 1",
        "--mpz-length-um 50",
        "--mpz-n-bins 80",
        "--da-phys 5e-6",
    )
    for token in required:
        assert token in text


def test_kernel_family_is_local_and_sha_pinned():
    text = RUNNER.read_text()
    assert (
        "$PARAMETER_SOURCE_ROOT/signed_kernel/"
        "v10_2_14_active_only_campaign_family.json"
    ) in text
    assert (
        "a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde"
        in text
    )


def test_validator_and_vendor_use_generic_local_input_language():
    validator = VALIDATOR.read_text()
    vendor = VENDOR.read_text()
    assert "--parameter-source-root" in validator
    assert "external_PF_runtime_dependency" in validator
    assert "--kernel-family-source" in vendor
    assert "kernel_family_vendored_byte_for_byte" in vendor
