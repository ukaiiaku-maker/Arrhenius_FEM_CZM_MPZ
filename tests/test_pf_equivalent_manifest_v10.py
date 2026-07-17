from arrhenius_fracture.pf_equivalent_material_manifest import (
    PF_SOURCE,
    load_material_manifest,
)


def test_exact_pf_manifests_load_without_refit():
    expected = {
        "ceramic": ("ceramic_restart02_candidate00", 12.718703137662922),
        "weakT": ("weakT_restart00_candidate00", 2.4387841773917582),
        "DBTT": ("DBTT_restart01_candidate05", 14.008733780506578),
    }
    for material, (candidate, sites) in expected.items():
        manifest = load_material_manifest(material, parameter_source=PF_SOURCE)
        assert manifest.candidate_id == candidate
        assert manifest.source_sites_per_system == sites
        assert manifest.parameter_source == PF_SOURCE
